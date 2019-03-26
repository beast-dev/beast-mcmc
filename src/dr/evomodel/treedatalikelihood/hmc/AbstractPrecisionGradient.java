/*
 * AbstractPrecisionGradient.java
 *
 * Copyright (c) 2002-2018 Alexei Drummond, Andrew Rambaut and Marc Suchard
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.evomodel.treedatalikelihood.hmc;

import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.CachedMatrixInverse;
import dr.inference.model.CompoundSymmetricMatrix;
import dr.inference.model.Likelihood;
import dr.inference.model.MatrixParameterInterface;
import dr.math.matrixAlgebra.Vector;
import dr.xml.Reportable;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

/**
 * @author Paul Bastide
 * @author Marc A. Suchard
 */

public abstract class AbstractPrecisionGradient implements GradientWrtParameterProvider, Reportable {

    private final GradientWrtPrecisionProvider gradientWrtPrecisionProvider;
    final Likelihood likelihood;
    final CompoundSymmetricMatrix compoundSymmetricMatrix;
    private final int dim;
    Parametrization parametrization;

    private final MatrixParameterInterface precision;
    private final MatrixParameterInterface variance;

    AbstractPrecisionGradient(GradientWrtPrecisionProvider gradientWrtPrecisionProvider,
                              Likelihood likelihood,
                              MatrixParameterInterface parameter) {

        this.precision = parameter;

        if (parameter instanceof CachedMatrixInverse) {

            this.compoundSymmetricMatrix = (CompoundSymmetricMatrix)
                    ((CachedMatrixInverse) parameter).getBaseParameter();
            this.variance = compoundSymmetricMatrix;
            this.parametrization = Parametrization.AS_VARIANCE;

        } else if (parameter instanceof CompoundSymmetricMatrix) {

            this.compoundSymmetricMatrix = (CompoundSymmetricMatrix) parameter;
            this.variance = new CachedMatrixInverse("", parameter);
            this.parametrization = Parametrization.AS_PRECISION;

        } else {
            throw new IllegalArgumentException("Unimplemented type");
        }

        assert compoundSymmetricMatrix.asCorrelation()
                : "PrecisionGradient can only be applied to a CompoundSymmetricMatrix with off-diagonal as correlation.";

        this.gradientWrtPrecisionProvider = gradientWrtPrecisionProvider;
        this.likelihood = likelihood;
        this.dim = parameter.getColumnDimension();
    }

    enum Parametrization {
        AS_PRECISION {
            @Override
            public double[] chainRule(double[] x, double[] vecP, double[] vecV) {
                return x; // Do nothing
            }
        },
        AS_VARIANCE {
            @Override
            public double[] chainRule(double[] x, double[] vecP, double[] vecV) {
                MultivariateChainRule ruleI = new MultivariateChainRule.InverseGeneral(vecP);
                return ruleI.chainGradient(x);
            }
        },
        AS_VARIANCE_DIAGONAL {
            @Override
            public double[] chainRule(double[] x, double[] vecP, double[] vecV) {
                MultivariateChainRule ruleI = new MultivariateChainRule.Inverse(vecP, vecV);
                return ruleI.chainGradient(x);
            }
        };

        abstract double[] chainRule(double[] x, double[] vecP, double[] vecV);
    }

    @Override
    public Likelihood getLikelihood() {
        return likelihood;
    }

    int getDimensionCorrelation() {
        return dim * (dim - 1) / 2;
    }

    int getDimensionDiagonal() {
        return dim;
    }

    @Override
    public double[] getGradientLogDensity() {

        // parameters
        double[] vecV = flatten(variance.getParameterAsMatrix());
        double[] vecP = flatten(precision.getParameterAsMatrix());
//        double[] diagQ = compoundSymmetricMatrix.getDiagonal();
//        double[] vecC = flatten(compoundSymmetricMatrix.getCorrelationMatrix());

        // Gradient w.r.t. precision
        double[] gradient = gradientWrtPrecisionProvider.getGradientWrtPrecision(vecV);

        // Handle inverse
        gradient = parametrization.chainRule(gradient, vecP, vecV);

        if (CHECK_GRADIENT) {
            System.err.println("Analytic at: \n" + new Vector(compoundSymmetricMatrix.getOffDiagonalParameter().getParameterValues())
                    + " " + new Vector(compoundSymmetricMatrix.getDiagonal()));
        }

        gradient = getGradientParameter(gradient);

        if (CHECK_GRADIENT) {
            System.err.println(checkNumeric(gradient));
        }

        return gradient;
    }

    String getReportString(double[] analytic, double[] numeric) {

        return getClass().getCanonicalName() + "\n" +
                "analytic: " + new Vector(analytic) +
                "\n" +
                "numeric: " + new Vector(numeric) +
                "\n";
    }

    String getReportString(double[] analytic, double[] numeric, double[] numericTrans) {

        return getClass().getCanonicalName() + "\n" +
                "analytic: " + new Vector(analytic) +
                "\n" +
                "numeric (no Cholesky): " + new Vector(numeric) +
                "\n" +
                "numeric (with Cholesky): " + new Vector(numericTrans) +
                "\n";
    }

    abstract String checkNumeric(double[] analytic);

    @Override
    public String getReport() {
        return checkNumeric(getGradientLogDensity());
    }

    abstract double[] getGradientParameter(double[] gradient);

    // Gradient w.r.t. correlation
    double[] getGradientCorrelation(double[] gradient) {

        gradient = compoundSymmetricMatrix.updateGradientOffDiagonal(gradient);

        return gradient;
    }

    // Gradient w.r.t. diagonal
    double[] getGradientDiagonal(double[] gradient) {

        gradient = compoundSymmetricMatrix.updateGradientDiagonal(gradient);

        return gradient;
    }

    public static double[] flatten(double[][] matrix) {
        int dim = 0;
        for (double[] vector : matrix) {
            dim += vector.length;
        }

        double[] result = new double[dim];

        int offset = 0;
        for (double[] vector : matrix) {
            System.arraycopy(vector, 0, result, offset, vector.length);
            offset += vector.length;
        }

        return result;
    }


    private static final boolean CHECK_GRADIENT = false;

    interface MultivariateChainRule {

        double[] chainGradient(double[] lhs);

        class Chain implements MultivariateChainRule {

            private final MultivariateChainRule[] rules;

            Chain(MultivariateChainRule[] rules) {
                this.rules = rules;
            }

            @Override
            public double[] chainGradient(double[] gradient) {

                for (MultivariateChainRule rule : rules) {
                    gradient = rule.chainGradient(gradient);
                }
                return gradient;
            }
        }

        class Inverse implements MultivariateChainRule {

            private final double[] vecP;
            private final double[] vecV;
            private final int dim;

            Inverse(double[] vecP, double[] vecV) {
                this.vecP = vecP;
                this.vecV = vecV;
                this.dim = (int) Math.sqrt(vecP.length);
            }

            @Override
            public double[] chainGradient(double[] lhs) {

                assert lhs.length == dim * dim;

                double[] gradient = new double[dim * dim];

                for (int i = 0; i < dim * dim; ++i) {

                    if (vecV[i] == 0 || Double.isNaN(vecV[i])) {
                        throw new RuntimeException("0 or NaN value in variance. check start value or use smaller step size for hmc");
                    }
                    gradient[i] = -lhs[i] * vecP[i] / vecV[i];
                }

                return gradient;
            }
        }

        class InverseGeneral implements MultivariateChainRule {

            private final DenseMatrix64F vecP;
            private final DenseMatrix64F temp;
            private final int dim;


            InverseGeneral(double[] vecP) {
                this.dim = (int) Math.sqrt(vecP.length);
                this.vecP = DenseMatrix64F.wrap(dim, dim, vecP);
                this.temp = new DenseMatrix64F(dim, dim);
            }

            @Override
            public double[] chainGradient(double[] lhs) {

                assert lhs.length == dim * dim;

                DenseMatrix64F gradient = new DenseMatrix64F(dim, dim);

                DenseMatrix64F LHS = DenseMatrix64F.wrap(dim, dim, lhs);
                CommonOps.mult(vecP, LHS, temp);
                CommonOps.mult(-1, temp, vecP, gradient);

                return gradient.getData();
            }
        }
    }
}

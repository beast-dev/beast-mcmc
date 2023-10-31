
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

import dr.evomodel.treedatalikelihood.continuous.ContinuousTraitGradientForBranch;
import dr.inference.model.*;
import dr.math.matrixAlgebra.Vector;
import dr.xml.Reportable;

/**
 * @author Paul Bastide
 * @author Marc A. Suchard
 */

public abstract class AbstractPrecisionGradient extends AbstractDiffusionGradient implements Reportable {

    private final GradientWrtPrecisionProvider gradientWrtPrecisionProvider;
    //    final Likelihood likelihood;
    protected final CompoundSymmetricMatrix compoundSymmetricMatrix;
    private final int dim;
    private Parametrization parametrization;

    private final MatrixParameterInterface precision;
    private final MatrixParameterInterface variance;

    AbstractPrecisionGradient(GradientWrtPrecisionProvider gradientWrtPrecisionProvider,
                              Likelihood likelihood,
                              MatrixParameterInterface parameter,
                              double upperBound, double lowerBound) {

        super(likelihood, upperBound, lowerBound);

        this.precision = parameter;

        if (parameter instanceof CachedMatrixInverse) {

            this.compoundSymmetricMatrix = (CompoundSymmetricMatrix)
                    ((CachedMatrixInverse) parameter).getBaseParameter();
            this.variance = compoundSymmetricMatrix;
            this.parametrization = Parametrization.AS_VARIANCE;

        } else if (parameter instanceof CompoundSymmetricMatrix) {

            this.compoundSymmetricMatrix = (CompoundSymmetricMatrix) parameter;
            this.variance = new CachedMatrixInverse("", precision);
            this.parametrization = Parametrization.AS_PRECISION;

        } else {
            throw new IllegalArgumentException("Unimplemented type");
        }

        assert compoundSymmetricMatrix.asCorrelation()
                : "PrecisionGradient can only be applied to a CompoundSymmetricMatrix with off-diagonal as correlation.";

        this.gradientWrtPrecisionProvider = gradientWrtPrecisionProvider;
//        this.likelihood = likelihood;
        this.dim = parameter.getColumnDimension();

    }

    enum Parametrization {
        AS_PRECISION {
//            @Override
//            public double[] chainRule(double[] x, double[] vecP, double[] vecV) {
//                return x; // Do nothing
//            }

            @Override
            public double[] getGradientWrtParameter(double[] gradient, double[] vecP, double[] vecV, GradientWrtPrecisionProvider gradientWrtPrecisionProvider) {
                return gradientWrtPrecisionProvider.getGradientWrtPrecision(vecV, gradient);
            }

            @Override
            void updateParameters(MatrixParameterInterface variance) {
                ((CachedMatrixInverse) variance).forceComputeInverse(); // ensure that variance is up to date
            }
        },
        AS_VARIANCE {
//            @Override
//            public double[] chainRule(double[] x, double[] vecP, double[] vecV) {
//                MultivariateChainRule ruleI = new MultivariateChainRule.InverseGeneral(vecP);
//                return ruleI.chainGradient(x);
//            }

            @Override
            public double[] getGradientWrtParameter(double[] gradient, double[] vecP, double[] vecV, GradientWrtPrecisionProvider gradientWrtPrecisionProvider) {
                return gradientWrtPrecisionProvider.getGradientWrtVariance(vecP, vecV, gradient);
            }

            @Override
            void updateParameters(MatrixParameterInterface variance) {
                // Do nothing
            }
//        },
//        AS_VARIANCE_DIAGONAL {
//            @Override
//            public double[] chainRule(double[] x, double[] vecP, double[] vecV) {
//                MultivariateChainRule ruleI = new MultivariateChainRule.Inverse(vecP, vecV);
//                return ruleI.chainGradient(x);
//            }
        };

//        abstract double[] chainRule(double[] x, double[] vecP, double[] vecV);

        abstract double[] getGradientWrtParameter(double[] gradient, double[] vecP, double[] vecV, GradientWrtPrecisionProvider gradientWrtPrecisionProvider);

        abstract void updateParameters(MatrixParameterInterface variance);
    }

//    @Override
//    public Likelihood getLikelihood() {
//        return likelihood;
//    }

    @Override
    public Parameter getRawParameter() {
        return precision;
    }

    @Override
    public ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient.DerivationParameter getDerivationParameter() {
        return ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient.DerivationParameter.WRT_VARIANCE;
    }

    protected int getDimensionCorrelation() {
        return dim * (dim - 1) / 2;
    }

    protected int getDimensionDiagonal() {
        return dim;
    }

    protected int getDimensionFull() {
        return dim * dim;
    }

    @Override
    public double[] getGradientLogDensity() {
        double[] gradient = (gradientWrtPrecisionProvider.getBranchSpecificGradient() == null) ? null : gradientWrtPrecisionProvider.getBranchSpecificGradient().getGradientLogDensity(); // Get gradient wrt variance
        return getGradientLogDensity(gradient);
    }

    public double[] getGradientLogDensity(double[] grad) {

        double[] gradient = new double[dim * dim];
        if (grad != null) System.arraycopy(grad, offset, gradient, 0, dim * dim);

        // parameters
        parametrization.updateParameters(variance);

        double[] vecV = flatten(variance.getParameterAsMatrix());
        double[] vecP = flatten(precision.getParameterAsMatrix());
        if (DEBUG) {
            System.err.println("vecV: " + new dr.math.matrixAlgebra.Vector(vecV));
            System.err.println("vecP: " + new dr.math.matrixAlgebra.Vector(vecP));
        }

        // Gradient w.r.t. precision
//        gradient = gradientWrtPrecisionProvider.getGradientWrtPrecision(vecV, gradient);
        gradient = parametrization.getGradientWrtParameter(gradient, vecP, vecV, gradientWrtPrecisionProvider);

        if (DEBUG) {
            System.err.println("Gradient Precision: " + new dr.math.matrixAlgebra.Vector(gradient));
        }

        // Handle inverse
//        gradient = parametrization.chainRule(gradient, vecP, vecV);

        if (CHECK_GRADIENT) {
            System.err.println("Analytic at: \n" + new Vector(compoundSymmetricMatrix.getOffDiagonalParameter().getParameterValues())
                    + " " + new Vector(compoundSymmetricMatrix.getDiagonal()));
        }

        gradient = getGradientParameter(gradient);

        if (CHECK_GRADIENT) {
            System.err.println(checkNumeric(gradient));
        }

        if (DEBUG) {
            System.err.println("Gradient Parameter: " + new dr.math.matrixAlgebra.Vector(gradient));
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

//    abstract String checkNumeric(double[] analytic);

//    @Override
//    public String getReport() {
//        return GradientWrtParameterProvider.getReportAndCheckForError(this,
//                Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY, TOLERANCE);
//    }

    abstract double[] getGradientParameter(double[] gradient);

    // Gradient w.r.t. correlation
    double[] getGradientCorrelation(double[] gradient) {

        return compoundSymmetricMatrix.updateGradientOffDiagonal(gradient);

    }

    // Gradient w.r.t. diagonal
    double[] getGradientDiagonal(double[] gradient) {

        return compoundSymmetricMatrix.updateGradientDiagonal(gradient);

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

    private static final boolean DEBUG = false;
}

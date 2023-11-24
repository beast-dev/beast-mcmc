/*
 * GaussianMarkovRandomField.java
 *
 * Copyright (c) 2002-2023 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.math.distributions.gp;

import dr.inference.distribution.RandomField;
import dr.inference.model.*;
import dr.math.distributions.RandomFieldDistribution;
import org.ejml.alg.dense.decomposition.chol.CholeskyDecompositionCommon_D64;
import org.ejml.data.DenseMatrix64F;
import org.ejml.factory.LinearSolverFactory;
import org.ejml.interfaces.linsol.LinearSolver;

import java.util.Arrays;
import java.util.List;

import static dr.math.distributions.MultivariateNormalDistribution.gradLogPdf;
import static dr.math.matrixAlgebra.missingData.MissingOps.invertAndGetDeterminant;

/**
 * @author Marc Suchard
 * @author Filippo Monti
 * //
 * Duvenaud DK, Nickisch H, Rasmussen C. Additive Gaussian processes. In Shawe-Taylor J, Zemel R, Bartlett P, Pereira F, Weinberger KQ (eds.), Advances in Neural Information Processing Systems, volume 24. Curran Associates, Inc., 2011.
 * URL <a href="https://proceedings.neurips.cc/paper/2011/file/4c5bde74a8f110656874902f07378009-Paper.pdf"/>
 */
public class AdditiveGaussianProcessDistribution extends RandomFieldDistribution {

    public static final String TYPE = "GaussianProcess";

    private final int order;

    private final int dim;
    private final Parameter meanParameter;
    private final Parameter nuggetParameter;
    private final List<BasisDimension> bases;
    private final RandomField.WeightProvider weightProvider;

    private final double[] mean;
    private final double[] tmp;
    private final DenseMatrix64F gramian;
    private final DenseMatrix64F precision;
    private final DenseMatrix64F variance;

    private double logDeterminant;
    
    private boolean meanKnown;
    private boolean precisionAndDeterminantKnown;
    private boolean gramianAndVarianceKnown;

    private static final boolean USE_CHOLESKY = true;

    public AdditiveGaussianProcessDistribution(String name,
                                               int order,
                                               int dim,
                                               Parameter meanParameter,
                                               Parameter nuggetParameter,
                                               List<BasisDimension> bases,
                                               RandomField.WeightProvider weightProvider) {
        super(name);

        if (order != 1) {
            throw new RuntimeException("Not yet implemented");
        }

        this.order = order;
        this.dim = dim;
        this.meanParameter = meanParameter;
        this.nuggetParameter = nuggetParameter;
        this.bases = bases;
        this.weightProvider = weightProvider;

        this.mean = new double[dim];
        this.tmp = new double[dim];
        this.gramian = new DenseMatrix64F(dim, dim);
        this.precision = new DenseMatrix64F(dim, dim);
        this.variance = new DenseMatrix64F(dim, dim);

        if (meanParameter != null) {
            addVariable(meanParameter);
        }

        if (nuggetParameter != null) {
            addVariable(nuggetParameter);
        }

        for (BasisDimension basis : bases) {
            GaussianProcessKernel kernel = basis.getKernel();
            if (kernel instanceof AbstractModel) {
                addModel((AbstractModel) kernel);
            }

            addVariable(basis.getDesignMatrix1());
            addVariable(basis.getDesignMatrix2());
        }

        if (weightProvider != null) {
            addModel(weightProvider);
        }
    }

    public int getOrder() { return order; }

    List<BasisDimension> getBases() { return bases; }

    private void computeGramianAndVariance() {

        computeAdditiveGramian(gramian, bases, order);

        // Add variance nugget as needed
        variance.set(gramian);

        if (nuggetParameter != null) {
            for (int i = 0; i < dim; ++i) {
                variance.add(i, i, getNugget(i));
            }
        }
    }

    private void computePrecisionAndDeterminant() {
        DenseMatrix64F variance = getVariance();
        if (USE_CHOLESKY) {
            LinearSolver<DenseMatrix64F> solver = LinearSolverFactory.symmPosDef(dim);
            if (!solver.setA(variance)) {
                throw new RuntimeException("Unable to decompose matrix");
            }

            solver.invert(precision);
            logDeterminant = 2 * computeLogDeterminantFromTriangularMatrix(
                    ((CholeskyDecompositionCommon_D64) solver.getDecomposition()).getT());
        } else {
            logDeterminant = invertAndGetDeterminant(variance, precision, true);
        }
    }

    private static double computeLogDeterminantFromTriangularMatrix(DenseMatrix64F T) {

        final int n = T.numCols;
        double[] t = T.getData();

        double sum = 0.0;
        int total = n * n;

        for(int i = 0; i < total; i += n + 1) {
            sum += Math.log(t[i]);
        }

        return sum;
    }

    private double[] getPrecision() {
        if (!precisionAndDeterminantKnown) {
            computePrecisionAndDeterminant();
            precisionAndDeterminantKnown = true;
        }
        return precision.getData();
    }

    private double getLogDeterminant() {
        if (!precisionAndDeterminantKnown) {
            computePrecisionAndDeterminant();
            precisionAndDeterminantKnown = true;
        }
        return logDeterminant;
    }

    @SuppressWarnings("unused")
    private DenseMatrix64F getGramian() {
        if (!gramianAndVarianceKnown) {
            computeGramianAndVariance();
            gramianAndVarianceKnown = true;
        }
        return gramian;
    }

    private DenseMatrix64F getVariance() {
        if (!gramianAndVarianceKnown) {
             computeGramianAndVariance();
             gramianAndVarianceKnown = true;
        }
        return variance;
    }

    private double getNugget(int i) {
        return nuggetParameter.getDimension() == 1 ?
                nuggetParameter.getParameterValue(0) :
                nuggetParameter.getParameterValue(i);
    }

    @Override
    public double[] getMean() {
        if (!meanKnown) {
            if (meanParameter == null) {
                Arrays.fill(mean, 0.0);
            } else if (meanParameter.getDimension() == 1) {
                Arrays.fill(mean, meanParameter.getParameterValue(0));
            } else {
                for (int i = 0; i < mean.length; ++i) {
                    mean[i] = meanParameter.getParameterValue(i);
                }
            }
            meanKnown = true;
        }
        return mean;
    }

    @Override
    public GradientProvider getGradientWrt(Parameter parameter) {
        throw new RuntimeException("Unknown parameter");
    }

    @Override
    public String getType() {
        return TYPE;
    }

    @Override
    public double[][] getScaleMatrix() {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public Variable<Double> getLocationVariable() {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public double logPdf(double[] x) {

        final double[] mean = getMean();
        final double[] diff = tmp;
        final double[] precision = getPrecision();

        for (int i = 0; i < dim; ++i) {
            diff[i] = x[i] - mean[i];
        }

        double exponent = 0.0;
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                exponent += diff[i] * precision[i * dim + j] * diff[j];
            }
        }
        
        return -0.5 * (dim * Math.log(2 * Math.PI) - getLogDeterminant()) - 0.5 * exponent;
    }

    @Override
    public int getDimension() { return dim; }

    @Override
    public double[] getGradientLogDensity(Object x) {
        return gradLogPdf((double[]) x, getMean(), getPrecision());
    }

    @Override
    public double[] getDiagonalHessianLogDensity(Object x) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public double[][] getHessianLogDensity(Object x) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public double[] nextRandom() {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        throw new IllegalArgumentException("Unknown model");
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    protected void storeState() {

    }

    @Override
    protected void restoreState() {

    }

    @Override
    protected void acceptState() { }

    public static class BasisDimension {

        private final GaussianProcessKernel kernel;
        private final DesignMatrix design1;
        private final DesignMatrix design2;

        public BasisDimension(GaussianProcessKernel kernel, DesignMatrix design1, DesignMatrix design2) {
            this.kernel = kernel;
            this.design1 = design1;
            this.design2 = design2;
        }

        public BasisDimension(GaussianProcessKernel kernel, DesignMatrix design) {
            this(kernel, design, design);
        }

        GaussianProcessKernel getKernel() { return kernel; }

        DesignMatrix getDesignMatrix1() { return design1; }

        DesignMatrix getDesignMatrix2() { return design2; }
    }

    public static void computeAdditiveGramian(DenseMatrix64F gramian,
                                              List<BasisDimension> bases,
                                              int order) {
        gramian.zero();

        final int rowDim = gramian.getNumRows();
        final int colDim = gramian.getNumCols();

        // 1st order contribution
        for (BasisDimension basis : bases) {
            final GaussianProcessKernel kernel = basis.getKernel();
            final DesignMatrix design1 = basis.getDesignMatrix1();
            final DesignMatrix design2 = basis.getDesignMatrix2();
            final double scale = kernel.getScale(); // TODO is this term necessary?  or scale only needed at the order-level

            for (int i = 0; i < rowDim; ++i) {
                for (int j = 0; j < colDim; ++j) {
                    double xi = design1.getParameterValue(i, 0); // TODO make generic dimension
                    double xj = design2.getParameterValue(j, 0); // TODO make generic dimension
                    gramian.add(i, j, scale * kernel.getCorrelation(xi, xj));
                }
            }
        }

        for (int n = 1; n < order; ++n) {
            // TODO higher-order terms via Newton-Girard formula
        }
    }
}

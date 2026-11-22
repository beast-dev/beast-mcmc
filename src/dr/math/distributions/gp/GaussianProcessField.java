/*
 * AdditiveGaussianProcessDistribution.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.math.distributions.gp;

import dr.inference.model.*;
import dr.math.distributions.RandomFieldDistribution;
import dr.xml.Reportable;
import org.ejml.alg.dense.decomposition.chol.CholeskyDecompositionCommon_D64;
import org.ejml.data.DenseMatrix64F;
import org.ejml.factory.LinearSolverFactory;
import org.ejml.interfaces.linsol.LinearSolver;

import java.util.Arrays;
import java.util.List;

import static dr.math.matrixAlgebra.missingData.MissingOps.invertAndGetDeterminant;

/**
 * @author Marc Suchard
 * @author Filippo Monti
 * //
 * Duvenaud DK, Nickisch H, Rasmussen C. Additive Gaussian processes. In Shawe-Taylor J, Zemel R, Bartlett P, Pereira F, Weinberger KQ (eds.), Advances in Neural Information Processing Systems, volume 24. Curran Associates, Inc., 2011.
 * URL <a href="https://proceedings.neurips.cc/paper/2011/file/4c5bde74a8f110656874902f07378009-Paper.pdf"/>
 */
public class GaussianProcessField extends RandomFieldDistribution implements Reportable {

    public static final String TYPE = "GaussianProcess";

    private final int order;

    private final int dim;
    private final Parameter orderVariance; // TODO should we reparametrize in terms of precision (for priors)?
    private final Parameter meanParameter;
    private final Parameter nuggetParameter;
    private final List<BasisDimension> bases;

    private double[] mean;
    private final double[] tmp;

    private double[] storedPrecision;
    private double[] storedMean;
    private double storedLogDeterminant = 0.0;
    private final DenseMatrix64F gramian;
    private final DenseMatrix64F precision;
    private final DenseMatrix64F variance;
    private final double[] tmpMatrix;

    private double[] diff;
    private double[] precisionDiff;

    private double logDeterminant;

    private boolean meanKnown = false;
    private boolean precisionAndDeterminantKnown = false;
    private boolean factorizationKnown = false;
    private boolean logDeterminantKnown = false;
    private boolean gramianAndVarianceKnown = false;
    private boolean precisionDiffKnown = false;

    private boolean storedMeanKnown = false;
    private boolean storedPrecisionAndDeterminantKnown = false;
    private boolean storedFactorizationKnown = false;
    private boolean storedLogDeterminantKnown = false;
    private boolean storedGramianAndVarianceKnown = false;
    private boolean storedPrecisionDiffKnown = false;
    private boolean storedFieldUpdated = true;
    private boolean covarianceChangedSinceStore = false;

    private LinearSolver<DenseMatrix64F> solver;

    private static final boolean USE_CHOLESKY = true;

    private Parameter field = null;
    private boolean fieldUpdated = true;

    private final boolean DEBUG = false;

    public void passParameter(Parameter field) {
        this.field = field;
        if (field != null) {
            addVariable(field);
        }
    }

    public GaussianProcessField(String name,
                                int dim,
                                Parameter orderVariance,
                                Parameter meanParameter,
                                Parameter nuggetParameter,
                                List<BasisDimension> bases) {
        super(name);

        this.order = orderVariance.getDimension();

        if (order != 1) {
            throw new RuntimeException("Not yet implemented");
        }

        this.dim = dim;
        this.orderVariance = orderVariance;
        this.meanParameter = meanParameter;
        this.nuggetParameter = nuggetParameter;
        this.bases = bases;

        this.mean = new double[dim]; // mu
        this.diff = new double[dim]; // x - mu
        this.gramian = new DenseMatrix64F(dim, dim);
        this.variance = new DenseMatrix64F(dim, dim); // K
        this.precision = new DenseMatrix64F(dim, dim); // P = K^{-1}
        this.storedPrecision = new double[dim * dim];
        this.storedMean = new double[dim]; // mu

        this.precisionDiff = new double[dim]; // -P (x - mu)
        this.storedDiff = new double[dim];
        this.storedPrecisionDiff = new double[dim];

        this.tmpMatrix = new double[dim * dim];
        this.tmp = new double[dim];

        this.solver = LinearSolverFactory.symmPosDef(dim);

        addVariable(orderVariance);

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
    }

    public int getOrder() {
        return order;
    }

    public Parameter getOrderVariance() {
        return orderVariance;
    }

    Parameter getFieldParameter() {
        return field;
    }

    Parameter getMeanParameter() {
        return meanParameter;
    }

    Parameter getNuggetParameter() {
        return nuggetParameter;
    }

    public List<BasisDimension> getBases() {
        return bases;
    }

    private void computeGramianAndVariance() {
        if (DEBUG) System.out.println("Computing Gramian and variance");

        computeAdditiveGramian(gramian, bases, orderVariance);

        // Add variance nugget as needed
        variance.set(gramian);

        if (nuggetParameter != null) {
            for (int i = 0; i < dim; ++i) {
                variance.add(i, i, getNugget(i));
            }
        }
    }

    private void computePrecisionAndDeterminant() {
        if (USE_CHOLESKY) {
            ensureFactorizationAndDeterminant();
            solver.invert(precision);
        } else {
            logDeterminant = invertAndGetDeterminant(getVariance(), precision, true);
        }
        precisionAndDeterminantKnown = true;
        logDeterminantKnown = true;
    }

    private void ensureFactorizationAndDeterminant() {
        if (factorizationKnown) {
            return;
        }

        if (!USE_CHOLESKY) {
            computePrecisionAndDeterminant();
            return;
        }

        if (!solver.setA(getVariance())) {
            if (DEBUG) System.out.println(bases.get(0).getKernel().getParameters());
            throw new RuntimeException("Unable to decompose matrix");
        }

        logDeterminant = 2 * computeLogDeterminantFromTriangularMatrix(
                ((CholeskyDecompositionCommon_D64) solver.getDecomposition()).getT());
        factorizationKnown = true;
        logDeterminantKnown = true;
    }

    private static double computeLogDeterminantFromTriangularMatrix(DenseMatrix64F T) {

        final int n = T.numCols;
        double[] t = T.getData();

        double sum = 0.0;
        int total = n * n;

        for (int i = 0; i < total; i += n + 1) {
            sum += Math.log(t[i]);
        }

        return sum;
    }

    protected double[] getPrecision() {
        return getPrecisionAsMatrix().getData();
    }

    protected DenseMatrix64F getPrecisionAsMatrix() {
        if (!precisionAndDeterminantKnown) {
            computePrecisionAndDeterminant();
        }
        return precision;
    }

    private double getLogDeterminant() {
        if (!logDeterminantKnown) {
            ensureFactorizationAndDeterminant();
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

    protected double getNugget(int i) {
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

    protected double[] getPrecisionDiff(double[] x) {
        ensurePrecisionDiff(x);
        return precisionDiff;
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
    public int getDimension() {
        return dim;
    }

    @Override
    public double logPdf(double[] x) {
        if (DEBUG) System.out.println("LogPdf: " + precisionAndDeterminantKnown + " " + gramianAndVarianceKnown);
        ensurePrecisionDiff(x);

        double exponent = 0.0;
        for (int i = 0; i < dim; ++i) {
            exponent -= diff[i] * precisionDiff[i];
        }

        if (field != null) {
            if (DEBUG) System.out.println(bases.get(0).getKernel().getParameters());
            if (DEBUG)
                System.out.println("First 4 values of precision: " + Arrays.toString(Arrays.copyOfRange(precision.getData(), 0, Math.min(4, diff.length))));
            if (DEBUG)
                System.out.println("First 4 values of diff: " + Arrays.toString(Arrays.copyOfRange(diff, 0, Math.min(4, diff.length))));
            if (DEBUG)
                System.out.println("First 4 values of precisionDiff: " + Arrays.toString(Arrays.copyOfRange(precisionDiff, 0, Math.min(4, precisionDiff.length))));
            if (DEBUG) System.out.println("exponent = " + exponent);
            if (DEBUG) System.out.println("logdet = " + getLogDeterminant());
        }
          return -0.5 * (dim * Math.log(2 * Math.PI) + getLogDeterminant()) - 0.5 * exponent;
    }

    @Override
    public double[] getGradientLogDensity(Object x) {
        if(DEBUG) System.out.println("GradientLogDensity");
        ensurePrecisionDiff((double[]) x);
        if (DEBUG) System.out.println("First 4 values of precisionDiff: " + Arrays.toString(Arrays.copyOfRange(precisionDiff, 0, Math.min(4, precisionDiff.length))));
        return precisionDiff;
    }

    public void ensurePrecisionDiff(double[] x) {
        if (field == null) fieldUpdated = true;

        if (fieldUpdated || !meanKnown) {
            if(DEBUG) System.out.println("Field updated");
            computeDiff(x);
            precisionDiffKnown = false;
            fieldUpdated = false;
        }

        if (!precisionDiffKnown) {
            computePrecisionDiff();
            precisionDiffKnown = true;
        }
    }

    public void computeDiff(double[] x) { // x - mu
        final double[] mean = getMean();
        for (int i = 0; i < dim; ++i) {
            diff[i] = x[i] - mean[i];
        }
    }

    public void computePrecisionDiff() { // -P (x - mu)
        final double[] precision = getPrecision();
        for (int i = 0; i < dim; ++i) {
            precisionDiff[i] = 0.0;
            for (int j = 0; j < dim; ++j) {
                precisionDiff[i] -= precision[i * dim + j] * diff[j];
            }
        }
    }

    private boolean containsKernel(Model model) {
        for (BasisDimension basis : bases) {
            if (model == basis.getKernel()) {
                return true;
            }
        }
        return false;
    }

    private void invalidateCovarianceCaches() {
        precisionAndDeterminantKnown = false;
        factorizationKnown = false;
        logDeterminantKnown = false;
        gramianAndVarianceKnown = false;
        precisionDiffKnown = false;
        covarianceChangedSinceStore = true;
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (containsKernel(model)) {
            if(DEBUG) System.out.println("Model changed event: kernel");
            invalidateCovarianceCaches();
            fireModelChanged();
        } else {
            throw new IllegalArgumentException("Unknown model");
        }
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
//        if(DEBUG) System.out.println("Variable changed event");
        if (variable == meanParameter) {
            meanKnown = false;
            fieldUpdated = true;
            precisionDiffKnown = false;
            fireModelChanged();
        } else if (variable == nuggetParameter) {
            invalidateCovarianceCaches();
            fireModelChanged();
        } else if (variable == orderVariance) {
            invalidateCovarianceCaches();
            fireModelChanged();
        } else if (variable == field) {
            if(DEBUG) System.out.println("Field changed event"); //TODO this is called for every entry wit compound parameter inside HMC
            fieldUpdated = true;
            precisionDiffKnown = false;
        } else {
            for (BasisDimension basis : bases) {
                if (variable == basis.getDesignMatrix1() || variable == basis.getDesignMatrix2()) {
                    invalidateCovarianceCaches();
                    fireModelChanged();
                    return;
                }
            }
        }
    }

    double[] storedPrecisionDiff;
    double[] storedDiff;

    @Override
    protected void storeState() {
        if(DEBUG) System.out.println("Storing states");
        if(DEBUG) System.out.println("Storing precision and determinant");

        storedMeanKnown = meanKnown;
        storedPrecisionAndDeterminantKnown = precisionAndDeterminantKnown;
        storedFactorizationKnown = factorizationKnown;
        storedLogDeterminantKnown = logDeterminantKnown;
        storedGramianAndVarianceKnown = gramianAndVarianceKnown;
        storedPrecisionDiffKnown = precisionDiffKnown;
        storedFieldUpdated = fieldUpdated;
        covarianceChangedSinceStore = false;

        if (precisionAndDeterminantKnown) {
            System.arraycopy(precision.getData(), 0, storedPrecision, 0, dim * dim);
        }
        storedLogDeterminant = logDeterminant;
        if (meanKnown) {
            System.arraycopy(mean, 0, storedMean, 0, dim);
        }
        if (field != null) {
            System.arraycopy(diff, 0, storedDiff, 0, dim);
            if (precisionDiffKnown) {
                System.arraycopy(precisionDiff, 0, storedPrecisionDiff, 0, dim);
            }
        }
    }

    @Override
    protected void restoreState() {
        if(DEBUG) System.out.println("Rejected state");
        final boolean covarianceChanged = covarianceChangedSinceStore;
        double[] temp;

        if (storedPrecisionAndDeterminantKnown) {
            temp = storedPrecision;
            storedPrecision = precision.getData();
            precision.setData(temp);
        }

        if (storedLogDeterminantKnown) {
            logDeterminant = storedLogDeterminant;
        }

        if (storedMeanKnown) {
            temp = storedMean;
            storedMean = mean;
            mean = temp;
        }

        if (field !=  null) {
            temp = storedDiff;
            storedDiff = diff;
            diff = temp;

            if (storedPrecisionDiffKnown) {
                temp = storedPrecisionDiff;
                storedPrecisionDiff = precisionDiff;
                precisionDiff = temp;
            }
        }

        meanKnown = storedMeanKnown;
        if (covarianceChanged) {
            gramianAndVarianceKnown = false;
            factorizationKnown = false;
            precisionAndDeterminantKnown = storedPrecisionAndDeterminantKnown;
            logDeterminantKnown = storedLogDeterminantKnown;
        } else {
            gramianAndVarianceKnown = storedGramianAndVarianceKnown || gramianAndVarianceKnown;
            factorizationKnown = storedFactorizationKnown || factorizationKnown;
            precisionAndDeterminantKnown = storedPrecisionAndDeterminantKnown || precisionAndDeterminantKnown;
            logDeterminantKnown = storedLogDeterminantKnown || logDeterminantKnown;
        }
        precisionDiffKnown = field != null && storedPrecisionDiffKnown;
        fieldUpdated = storedFieldUpdated;
        covarianceChangedSinceStore = false;

        if (DEBUG) System.out.println("Restored states");
        if (DEBUG) System.out.println(bases.get(0).getKernel().getParameters());
        if (DEBUG)
            System.out.println("First 4 values of precision: " + Arrays.toString(Arrays.copyOfRange(precision.getData(), 0, Math.min(4, diff.length))));
        if (DEBUG)
            System.out.println("First 4 values of diff: " + Arrays.toString(Arrays.copyOfRange(diff, 0, Math.min(4, diff.length))));
        if (DEBUG)
            System.out.println("First 4 values of precisionDiff: " + Arrays.toString(Arrays.copyOfRange(precisionDiff, 0, Math.min(4, precisionDiff.length))));
        if (DEBUG) System.out.println("logdet = " + getLogDeterminant());
    }

    @Override
    protected void acceptState() {
        if(DEBUG) System.out.println("Accepting state");
        covarianceChangedSinceStore = false;
    }

    @Override
    public GradientProvider getGradientWrt(Parameter parameter) {
        return GaussianProcessFieldGradient.createProvider(this, parameter);
    }

    @Override
    public double[] getDiagonalHessianLogDensity(Object x) {
        return GaussianProcessFieldHessian.getDiagonalHessianLogDensity(this);
    }

    @Override
    public double[][] getHessianLogDensity(Object x) {
        return GaussianProcessFieldHessian.getHessianLogDensity(this);
    }

    @Override
    public double[] nextRandom() {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public String getReport() {
        final double[] mean = getMean();
        final double[] precision = getPrecision();
        final double[] precisionDiff = getPrecisionDiff(new double[dim]);

        StringBuilder sb = new StringBuilder();
        sb.append("predictionDiff: [");
        for (double value : precisionDiff) {
            sb.append(" ").append(value);
        }
        sb.append(" ]\n");
        sb.append("mean: [");
        for (double value : mean) {
            sb.append(" ").append(value);
        }
        sb.append(" ]\n");
        sb.append("precision: [");
        for (double value : precision) {
            sb.append(" ").append(value);
        }
        sb.append(" ]\n");
        return sb.toString();
    }
    
    public static void computeAdditiveGramian(DenseMatrix64F gramian,
                                              List<BasisDimension> bases,
                                              Parameter orderVariance) {

        gramian.zero();
        // 1st order contribution
        for (BasisDimension basis : bases) {
            basis.addFirstOrderKernelComponent(gramian);
        }

        // TODO higher-order terms via Newton-Girard formula
        @SuppressWarnings("unused")
        final int order = orderVariance.getDimension();
//        for (int n = 1; n < order; ++n) { }
    }

}

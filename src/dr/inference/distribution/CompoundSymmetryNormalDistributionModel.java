package dr.inference.distribution;

import dr.inference.model.*;
import dr.inferencexml.distribution.MultivariateNormalDistributionModelParser;
import dr.math.distributions.GaussianProcessRandomGenerator;
import dr.math.distributions.MultivariateNormalDistribution;
import dr.math.matrixAlgebra.Matrix;

/**
 * A class that acts as a model for normally distributed data with zero-mean and compound symmetry covariance.
 *
 * @author Zhenyu Zhang
 */
public class CompoundSymmetryNormalDistributionModel extends AbstractModel implements
        ParametricMultivariateDistributionModel, GaussianProcessRandomGenerator, GradientProvider, HessianProvider {

    private final int dim;
    private final Parameter marginal;
    private final Parameter rho;

    private MultivariateNormalDistribution distribution;
    private MultivariateNormalDistribution storedDistribution;
    private double[][] storedPrecision;

    private boolean distributionKnown;
    private boolean storedDistributionKnown;

    public CompoundSymmetryNormalDistributionModel(int dim,
                                                   Parameter marginalParameter,
                                                   Parameter rho) {
        super(MultivariateNormalDistributionModelParser.NORMAL_DISTRIBUTION_MODEL);

        this.dim = dim;

        this.marginal = marginalParameter;
        addVariable(marginalParameter);

        this.rho = rho;
        addVariable(rho);

        storedPrecision = getPrecisionMatrix();
        distribution = createNewDistribution();
        distributionKnown = true;
    }

    private MultivariateNormalDistribution createNewDistribution() {
        return new MultivariateNormalDistribution(new double[getDimension()], storedPrecision);
    }

    @Override
    public double[] nextRandom() {
        checkDistribution();
        throw new RuntimeException("Not yet implemented");
    }

    private void checkDistribution() {
        if (!distributionKnown) {
            distribution = createNewDistribution();
            distributionKnown = true;
        }
    }

    @Override
    public double logPdf(Object x) {
        checkDistribution();
        return distribution.logPdf(x);
    }

    @Override
    public Variable<Double> getLocationVariable() {
        return null;
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {

    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {

    }

    @Override
    protected void storeState() {
        storedDistribution = distribution;
        storedDistributionKnown = distributionKnown;
    }

    protected void restoreState() {
        distributionKnown = storedDistributionKnown;
        distribution = storedDistribution;
    }

    @Override
    protected void acceptState() {
    }

    @Override
    public double[] getDiagonalHessianLogDensity(Object x) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public double[][] getHessianLogDensity(Object x) {
        throw new RuntimeException("Not yet implemented");
    }

    public double[] getPrecisionVectorProduct(double[] x) {
        checkDistribution();

        final double[] product = new double[dim];

        for (int i = 0; i < dim; ++i) {
            double sum = 0;
            for (int j = 0; j < dim; ++j) {
                sum += storedPrecision[i][j] * x[j];
            }
            product[i] = sum;
        }
        return product;
    }

    public double[] getDiagonal() {
        checkDistribution();
        double[] diagonal = new double[dim];
        for (int i = 0; i < dim; i++) {
            diagonal[i] = storedPrecision[i][i];
        }
        return diagonal;
    }

    @Override
    public double[] getGradientLogDensity(Object x) {
        return distribution.getGradientLogDensity(x);
    }

    @Override
    public Likelihood getLikelihood() {
        return null;
    }

    @Override
    public int getDimension() {
        return dim;
    }

    @Override
    public double[][] getPrecisionMatrix() {
        // with the formula precision = (1 / (1 - rho)) * [I - (rho /(1 - rho + rho * d)) * J], J is all one matrix.
        double[][] precision = new double[dim][dim];
        double rhoValue = rho.getParameterValue(0);
        double marginalValue = marginal.getParameterValue(0);

        if (marginalValue != 1) {
            throw new RuntimeException("not yet implemented!");
        }

        double denom = (1 - rhoValue) * (1 - rhoValue + rhoValue * dim);
        double diagonalTerm = (1 - 2 * rhoValue + dim * rhoValue) / denom;
        double offDiagonalTerm = (-rhoValue) / denom;

        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                if (i != j) {
                    precision[i][j] = offDiagonalTerm;
                } else {
                    precision[i][j] = diagonalTerm;
                }
            }
        }
        return precision;
    }

    @Override
    public double logPdf(double[] x) {
        checkDistribution();
        return distribution.logPdf(x);
    }

    @Override
    public double[][] getScaleMatrix() {
        throw new RuntimeException("Not yet implemented");
    }

    public double[] getPrecisionColumn(int index) {
        checkDistribution();
        double[] column = new double[dim];
        for (int i = 0; i < dim; i++) {
            column[i] = storedPrecision[i][index];
        }
        return column;
    }

    @Override
    public double[] getMean() {
        return new double[dim];
    }

    @Override
    public String getType() {
        return distribution.getType();
    }
}

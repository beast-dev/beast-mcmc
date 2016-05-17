package dr.inference.operators;

import dr.evomodel.continuous.FullyConjugateMultivariateTraitLikelihood;
import dr.inference.model.LatentFactorModel;
import dr.inference.model.MatrixParameterInterface;
import dr.math.MathUtils;
import dr.math.distributions.MultivariateNormalDistribution;
import dr.math.matrixAlgebra.Matrix;
import dr.math.matrixAlgebra.SymmetricMatrix;

/**
 * Created by max on 5/16/16.
 */
public class FactorTreeGibbsOperator extends SimpleMCMCOperator implements GibbsOperator{
    LatentFactorModel lfm;
    double pathParameter = 1;
    FullyConjugateMultivariateTraitLikelihood tree;
    MatrixParameterInterface factors;
    MatrixParameterInterface errorPrec;

    public FactorTreeGibbsOperator(double weight, LatentFactorModel lfm, FullyConjugateMultivariateTraitLikelihood tree){
        setWeight(weight);
        this.tree = tree;
        this.lfm = lfm;
        this.factors = lfm.getFactors();
        errorPrec = lfm.getColumnPrecision();
    }

    @Override
    public int getStepCount() {
        return 0;
    }

    @Override
    public String getPerformanceSuggestion() {
        return null;
    }

    @Override
    public String getOperatorName() {
        return null;
    }

    @Override
    public double doOperation() throws OperatorFailedException {
        int column = MathUtils.nextInt(factors.getColumnDimension());
        MultivariateNormalDistribution mvn = getMVN(column);
        double[] draw = (double[]) mvn.nextRandom();
        for (int i = 0; i < factors.getRowDimension(); i++) {
            factors.setParameterValueQuietly(i, column, draw[i]);
        }
        factors.fireParameterChangedEvent(column * factors.getRowDimension() , null);

        return 0;
    }

    MultivariateNormalDistribution getMVN(int column){
        double[][] precision = getPrecision(column);
        double[] mean = getMean(column, precision);
        return new MultivariateNormalDistribution(mean, precision);
    }

    double[][] getPrecision(int column){
        double [][] treePrec = tree.getConditionalPrecision(column);

        for (int i = 0; i < lfm.getLoadings().getColumnDimension(); i++) {
            for (int j = i; j < lfm.getLoadings().getColumnDimension(); j++) {
                for (int k = 0; k < lfm.getLoadings().getRowDimension(); k++) {
                    treePrec[i][j] += lfm.getLoadings().getParameterValue(k, i) * errorPrec.getParameterValue(k, k) * lfm.getLoadings().getParameterValue(k, j) * pathParameter;
                    treePrec[j][i] = treePrec[i][j];
                }
            }
        }
        return treePrec;
    }

    double[] getMean(int column, double[][] precision){
        Matrix variance = (new SymmetricMatrix(precision)).inverse();
        double[] midMean = new double[lfm.getLoadings().getRowDimension()];
        double[] condMean = tree.getConditionalMean(column);
        double[][] condPrec = tree.getConditionalPrecision(column);
        for (int i = 0; i < midMean.length; i++) {
            for (int j = 0; j < midMean.length; j++) {
                midMean [i] += condPrec[i][j] * condMean[j];
            }
        }
        for (int i = 0; i < lfm.getLoadings().getRowDimension(); i++) {
            for (int j = 0; j < lfm.getLoadings().getColumnDimension(); j++) {
                midMean[j] += lfm.getScaledData().getParameterValue(j, column) * errorPrec.getParameterValue(i,i) * lfm.getLoadings().getParameterValue(i, j);
            }
        }
        double[] mean = new double[midMean.length];
        for (int i = 0; i < mean.length; i++) {
            for (int j = 0; j < mean.length; j++) {
                mean[i] += variance.component(i, j) * midMean[j];
            }

        }
        return mean;
    }

    @Override
    public void setPathParameter(double beta) {
        pathParameter = beta;
    }
}

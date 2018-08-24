package dr.inference.operators;

import dr.evomodel.continuous.GibbsSampleFromTreeInterface;
import dr.inference.model.LatentFactorModel;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Parameter;
import dr.math.MathUtils;
import dr.math.distributions.MultivariateNormalDistribution;
import dr.math.matrixAlgebra.Matrix;
import dr.math.matrixAlgebra.SymmetricMatrix;

/**
 * Created by max on 5/16/16.
 */
public class FactorTreeGibbsOperator extends SimpleMCMCOperator implements PathDependentOperator, GibbsOperator {

    private final LatentFactorModel lfm;
    private double pathParameter = 1;
    private final GibbsSampleFromTreeInterface tree;
    private final GibbsSampleFromTreeInterface workingTree;
    private final MatrixParameterInterface factors;
    private final MatrixParameterInterface errorPrec;
    private final boolean randomScan;
    private final Parameter missingIndicator;

    public FactorTreeGibbsOperator(double weight, LatentFactorModel lfm, GibbsSampleFromTreeInterface tree, Boolean randomScan){
        setWeight(weight);
        this.tree = tree;
        this.lfm = lfm;
        this.factors = lfm.getFactors();
        errorPrec = lfm.getColumnPrecision();
        this.randomScan = randomScan;
        this.workingTree = null;
        missingIndicator = lfm.getMissingIndicator();
    }
    
    @Override
    public String getPerformanceSuggestion() {
        return null;
    }

    @Override
    public String getOperatorName() {
        return "Factor Tree Gibbs Operator";
    }

    @Override
    public double doOperation() {
        if(randomScan){
            int column = MathUtils.nextInt(factors.getColumnDimension());
            MultivariateNormalDistribution mvn = getMVN(column);
            double[] draw = (double[]) mvn.nextRandom();
            for (int i = 0; i < factors.getRowDimension(); i++) {
                factors.setParameterValue(i, column, draw[i]);
            }
        }
        else{
            for (int i = 0; i < factors.getColumnDimension(); i++) {
                MultivariateNormalDistribution mvn = getMVN(i);
                double[] draw = (double[]) mvn.nextRandom();
                for (int j = 0; j < factors.getRowDimension(); j++) {
                    factors.setParameterValue(j, i, draw[j]);
                }
            }
        }

        return 0;
    }

    MultivariateNormalDistribution getMVN(int column){
        double[][] precision = getPrecision(column);
        double[] mean = getMean(column, precision);
        return new MultivariateNormalDistribution(mean, precision);
    }

    double[][] getPrecision(int column){
        double [][] treePrec = getTreePrec(column);
        for (int i = 0; i < lfm.getLoadings().getColumnDimension(); i++) {
            for (int j = i; j < lfm.getLoadings().getColumnDimension(); j++) {
                for (int k = 0; k < lfm.getLoadings().getRowDimension(); k++) {
                    if(missingIndicator == null || missingIndicator.getParameterValue(column * lfm.getLoadings().getRowDimension() + k) != 1){
                    treePrec[i][j] += lfm.getLoadings().getParameterValue(k, i) * errorPrec.getParameterValue(k, k) * lfm.getLoadings().getParameterValue(k, j) * pathParameter;}
                }
                treePrec[j][i] = treePrec[i][j];
            }
        }
        return treePrec;
    }

    double[] getMean(int column, double[][] precision){
        Matrix variance = (new SymmetricMatrix(precision)).inverse();
        double[] midMean = new double[lfm.getLoadings().getColumnDimension()];
        double[] condMean = getTreeMean(column);
        double[][] condPrec = getTreePrec(column);
        for (int i = 0; i < midMean.length; i++) {
//            for (int j = 0; j < midMean.length; j++) {
                midMean [i] += condPrec[i][i] * condMean[i];
//            }
        }
        for (int i = 0; i < lfm.getLoadings().getRowDimension(); i++) {
            for (int j = 0; j < lfm.getLoadings().getColumnDimension(); j++) {
                if(missingIndicator == null || missingIndicator.getParameterValue(column * lfm.getScaledData().getRowDimension() + i) != 1)
                    midMean[j] += lfm.getScaledData().getParameterValue(i, column) * errorPrec.getParameterValue(i,i) * lfm.getLoadings().getParameterValue(i, j) * pathParameter;
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

    public double[][] getTreePrec(int column){
        double answerFactor = tree.getPrecisionFactor(column);
        double[][] answer = new double[factors.getRowDimension()][factors.getRowDimension()];

        for (int i = 0; i < factors.getRowDimension(); i++) {
            answer[i][i] = answerFactor;
        }

        if (workingTree != null) {
            double[][] temp = workingTree.getConditionalPrecision(column);
            for (int i = 0; i < answer.length; i++) {
                for (int j = 0; j < answer.length; j++) {
                    answer[i][j] = answer[i][j] * pathParameter + temp[i][j] * (1 - pathParameter);
                }

            }
        }

        return answer;
    }

    public double[] getTreeMean(int column){
        double[] answer = tree.getConditionalMean(column);

        if (workingTree != null) {
            double[] temp = workingTree.getConditionalMean(column);
            for (int i = 0; i < answer.length; i++) {
                answer[i] = answer[i] * pathParameter + temp[i] * (1 - pathParameter);
            }
        }

        return answer;
    }

    @Override
    public void setPathParameter(double beta) {
        pathParameter = beta;
    }
}

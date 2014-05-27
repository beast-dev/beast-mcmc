package dr.inference.operators;

import dr.evomodel.continuous.LatentFactorModel;
import dr.inference.distribution.DistributionLikelihood;
import dr.inference.model.MatrixParameter;
import dr.inference.model.Parameter;
import dr.math.distributions.NormalDistribution;
import dr.math.matrixAlgebra.Matrix;

/**
 * Created with IntelliJ IDEA.
 * User: max
 * Date: 5/23/14
 * Time: 2:23 PM
 * To change this template use File | Settings | File Templates.
 */
public class LoadingsGibbsOperator extends SimpleMCMCOperator implements GibbsOperator{
    NormalDistribution prior;
    LatentFactorModel LFM;

        public LoadingsGibbsOperator(LatentFactorModel LFM, DistributionLikelihood prior, double weight){
        setWeight(weight);

        this.prior=(NormalDistribution) prior.getDistribution();
        this.LFM=LFM;
    }

    private Matrix truncateMatrixParameter(MatrixParameter full, int i){
        double[][] answerArray=new double[i][full.getColumnDimension()];
        for (int j = 0; j <i ; j++) {
        answerArray[j]=full.getRowValues(j);
        }
//        System.out.println(answerArray.length);

        return new Matrix(answerArray);
    }

    @Override
    public int getStepCount() {
        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public String getPerformanceSuggestion() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public String getOperatorName() {
        return "loadingsGibbsOperator";  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public double doOperation() throws OperatorFailedException {
        Matrix tempFactors = null;
        int size = LFM.getLoadings().getColumnDimension();
        for (int i = 0; i < size; i++) {
            if (i < size) {
                tempFactors = truncateMatrixParameter(LFM.getFactors(), i+1);
//                System.out.println(tempFactors.rows());
//                System.out.println(tempFactors.columns());
//                System.out.println(tempFactors);
//                System.out.println(new Matrix(LFM.getFactors().getParameterAsMatrix()));
            }
        }


        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }
}

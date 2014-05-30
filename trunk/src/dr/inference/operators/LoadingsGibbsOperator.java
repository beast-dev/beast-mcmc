package dr.inference.operators;

import dr.evomodel.continuous.LatentFactorModel;
import dr.inference.distribution.DistributionLikelihood;
import dr.inference.model.MatrixParameter;
import dr.inference.model.Parameter;
import dr.math.distributions.MultivariateNormalDistribution;
import dr.math.distributions.NormalDistribution;
import dr.math.matrixAlgebra.IllegalDimension;
import dr.math.matrixAlgebra.Matrix;
import dr.math.matrixAlgebra.Vector;

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

    private Matrix getPrecision(int i){
        Matrix factors=null;
        int size=LFM.getFactorDimension();
        Matrix answer=null;
        if(i<size){
            factors=truncateMatrixParameter(LFM.getFactors(),i+1);
            try {
                answer= Matrix.buildIdentityTimesElementMatrix(i + 1, 1 / (prior.getSD() * prior.getSD())).add(factors.productWithTransposed(factors).product(LFM.getColumnPrecision().getParameterValue(i, i)));
            } catch (IllegalDimension illegalDimension) {
                illegalDimension.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            }
        }
        else{
                answer=new Matrix(LFM.getFactors().getParameterAsMatrix());
                try {
                    answer= Matrix.buildIdentityTimesElementMatrix(i, 1/(prior.getSD()*prior.getSD())).add(factors.productWithTransposed(factors).product(LFM.getColumnPrecision().getParameterValue(i, i)));
                } catch (IllegalDimension illegalDimension) {
                    illegalDimension.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
                }
            }
        return answer;
    }

    private Vector getMean(int i, Matrix precision){
        Matrix factors=null;
        int size=LFM.getFactorDimension();
        Vector answer=null;
        Vector dataColumn=null;
        Vector priorVector=null;
        Vector temp=null;
        if(i<size){
            factors=truncateMatrixParameter(LFM.getFactors(), i+1);
            dataColumn=new Vector(LFM.getScaledData().toComponents()[i]);
            priorVector=Vector.buildOneTimesElementVector(i+1, prior.getMean()/(prior.getSD()*prior.getSD()));
            try {
                answer=precision.inverse().product(priorVector.add(factors.product(dataColumn)));
            } catch (IllegalDimension illegalDimension) {
                illegalDimension.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            }
        }
        else{
            factors=new Matrix(LFM.getFactors().getParameterAsMatrix());
            dataColumn=new Vector(LFM.getScaledData().toComponents()[i]);
            priorVector=Vector.buildOneTimesElementVector(size, prior.getMean()/(prior.getSD()*prior.getSD()));
            try {
                answer=precision.inverse().product(priorVector.add(factors.product(dataColumn)));
            } catch (IllegalDimension illegalDimension) {
                illegalDimension.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            }
        }


        return answer;
    }

    private void copy(int i, double[] random){
        Parameter changing=LFM.getLoadings().getParameter(i);
        for (int j = 0; j <random.length ; j++) {
            changing.setParameterValueQuietly(j,random[j]);
        }
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
        double[] draws=null;
        int size = LFM.getLoadings().getColumnDimension();
        for (int i = 0; i < size; i++) {
            Matrix precision=getPrecision(i);
            Vector mean=getMean(i, precision);
            draws= MultivariateNormalDistribution.nextMultivariateNormalPrecision(mean.toComponents(),precision.toComponents());
            while(draws[i]<0){
                draws= MultivariateNormalDistribution.nextMultivariateNormalPrecision(mean.toComponents(),precision.toComponents());
            }
            copy(i, draws);
        }
        LFM.getLoadings().fireParameterChangedEvent();

        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }
}

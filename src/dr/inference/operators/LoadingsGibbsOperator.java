package dr.inference.operators;

import dr.evomodel.continuous.LatentFactorModel;
import dr.inference.distribution.DistributionLikelihood;
import dr.inference.model.MatrixParameter;
import dr.inference.model.Parameter;
import dr.math.distributions.MultivariateNormalDistribution;
import dr.math.distributions.NormalDistribution;
import dr.math.matrixAlgebra.SymmetricMatrix;

import java.util.ArrayList;
import java.util.ListIterator;

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
    ArrayList<double[][]> precisionArray;
    ArrayList<double[]> meanMidArray;
    ArrayList<double[]> meanArray;
    MatrixParameter[] vectorProductAnswer;
    MatrixParameter[] priorMeanVector;

    double priorPrecision;
    double priorMeanPrecision;

        public LoadingsGibbsOperator(LatentFactorModel LFM, DistributionLikelihood prior, double weight){
        setWeight(weight);

        this.prior=(NormalDistribution) prior.getDistribution();
        this.LFM=LFM;
        precisionArray=new ArrayList<double[][]>();
            double[][] temp;
        for (int i = 0; i < LFM.getFactorDimension() ; i++) {
            temp=new double[i+1][i+1];
            precisionArray.add(temp);
            }


            meanArray=new ArrayList<double[]>();
            meanMidArray=new ArrayList<double[]>();
            double[] tempMean;
            for (int i = 0; i < LFM.getFactorDimension() ; i++) {
                tempMean=new double[i+1];
                meanArray.add(tempMean);
            }

            for (int i = 0; i < LFM.getFactorDimension() ; i++) {
                tempMean=new double[i+1];
                meanMidArray.add(tempMean);
            }


            vectorProductAnswer=new MatrixParameter[LFM.getLoadings().getRowDimension()];
            for (int i = 0; i <vectorProductAnswer.length ; i++) {
                vectorProductAnswer[i]=new MatrixParameter(null);
                vectorProductAnswer[i].setDimensions(i+1, 1);
            }

        priorMeanVector=new MatrixParameter[LFM.getLoadings().getRowDimension()];
            for (int i = 0; i <priorMeanVector.length ; i++) {
                priorMeanVector[i]=new MatrixParameter(null, i+1, 1, this.prior.getMean()/(this.prior.getSD()*this.prior.getSD()));


            }
            priorPrecision= 1/(this.prior.getSD()*this.prior.getSD());
            priorMeanPrecision=this.prior.getMean()*priorPrecision;
    }

    private void getPrecisionOfTruncated(MatrixParameter full, int newRowDimension, double[][] answer){

//        MatrixParameter answer=new MatrixParameter(null);
//        answer.setDimensions(this.getRowDimension(), Right.getRowDimension());
//        System.out.println(answer.getRowDimension());
//        System.out.println(answer.getColumnDimension());

        int p = full.getColumnDimension();
        for (int i = 0; i < newRowDimension; i++) {
            for (int j = i; j < newRowDimension; j++) {
                double sum = 0;
                for (int k = 0; k < p; k++)
                    sum += full.getParameterValue(i, k) * full.getParameterValue(j,k);
                answer[i][j]=sum*LFM.getColumnPrecision().getParameterValue(newRowDimension,newRowDimension);
                if(i==j) {
                    answer[i][j] += priorPrecision;
                }
                else{
                    answer[j][i]=answer[i][j];
                }
                }
            }
        }


    private void getTruncatedMean(int newRowDimension, int dataColumn, double[][] variance, double[] midMean, double[] mean){

//        MatrixParameter answer=new MatrixParameter(null);
//        answer.setDimensions(this.getRowDimension(), Right.getRowDimension());
//        System.out.println(answer.getRowDimension());
//        System.out.println(answer.getColumnDimension());
        MatrixParameter data=LFM.getScaledData();
        MatrixParameter Left=LFM.getFactors();
        int p = Left.getColumnDimension();
        for (int i = 0; i < newRowDimension; i++) {
            double sum = 0;
                for (int k = 0; k < p; k++)
                    sum += Left.getParameterValue(i, k) * data.getParameterValue(dataColumn, k);
                sum=sum*LFM.getColumnPrecision().getParameterValue(i,i);
                sum+=priorMeanPrecision;
                midMean[i]=sum;
            }
        for (int i = 0; i < newRowDimension; i++) {
            double sum = 0;
            for (int k = 0; k < newRowDimension; k++)
                sum += variance[i][k] * midMean[k];
            mean[i]=sum;
        }

    }

    private void getPrecision(int i, double[][] answer){
        int size=LFM.getFactorDimension();
        if(i<size){
            getPrecisionOfTruncated(LFM.getFactors(), i + 1, answer);
        }
        else{
             getPrecisionOfTruncated(LFM.getFactors(), size, answer);
            }
    }

    private void getMean(int i, double[][] variance, double[] midMean, double[] mean){
//        Matrix factors=null;
        int size=LFM.getFactorDimension();
//        double[] scaledDataColumn=LFM.getScaledData().getRowValues(i);
//        Vector dataColumn=null;
//        Vector priorVector=null;
//        Vector temp=null;
//        Matrix data=new Matrix(LFM.getScaledData().getParameterAsMatrix());
        if(i<size){
            getTruncatedMean(i + 1, i, variance, midMean, mean);
//            dataColumn=new Vector(data.toComponents()[i]);
//            try {
//                answer=precision.inverse().product(new Matrix(priorMeanVector[i].add(vectorProductAnswer[i]).getParameterAsMatrix()));
//            } catch (IllegalDimension illegalDimension) {
//                illegalDimension.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
//            }
        }
        else{
            getTruncatedMean(size, i, variance, midMean, mean);
//            dataColumn=new Vector(data.toComponents()[i]);
//            try {
//                answer=precision.inverse().product(new Matrix(priorMeanVector[size-1].add(vectorProductAnswer[size-1]).getParameterAsMatrix()));
//            } catch (IllegalDimension illegalDimension) {
//                illegalDimension.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
//            }
        }

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
        ListIterator<double[][]> currentPrecision=precisionArray.listIterator();
        ListIterator<double[]> currentMidMean=meanMidArray.listIterator();
        ListIterator<double[]> currentMean=meanArray.listIterator();
        double[] draws;
        double[][] precision=null;
        double[][] variance;
        double[] midMean=null;
        double[] mean=null;
        int size = LFM.getLoadings().getColumnDimension();
        for (int i = 0; i < size; i++) {
            if(currentPrecision.hasNext()){
                precision=currentPrecision.next();}
            if(currentMidMean.hasNext())
            {midMean=currentMidMean.next();
            }
            if(currentMean.hasNext()){
                mean=currentMean.next();
            }
            getPrecision(i, precision);
            variance=(new SymmetricMatrix(precision)).inverse().toComponents();
            getMean(i, variance, mean, midMean);

            draws= MultivariateNormalDistribution.nextMultivariateNormalVariance(mean, variance);
            if(i<draws.length){
                while(draws[i]<0){
                    draws= MultivariateNormalDistribution.nextMultivariateNormalVariance(mean, variance);
                }
            }
            copy(i, draws);
        }
        LFM.getLoadings().fireParameterChangedEvent();

        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }
}

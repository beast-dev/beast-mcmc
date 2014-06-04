package dr.inference.operators;

import dr.evomodel.continuous.LatentFactorModel;
import dr.inference.distribution.DistributionLikelihood;
import dr.inference.model.DiagonalMatrix;
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
    DiagonalMatrix[] singleElement;
    MatrixParameter[] innerProductAnswer;
    MatrixParameter[] vectorProductAnswer;
    MatrixParameter[] priorMeanVector;

        public LoadingsGibbsOperator(LatentFactorModel LFM, DistributionLikelihood prior, double weight){
        setWeight(weight);

        this.prior=(NormalDistribution) prior.getDistribution();
        this.LFM=LFM;
        singleElement=new DiagonalMatrix[LFM.getLoadings().getRowDimension()];
            for (int i = 0; i <singleElement.length ; i++) {
                singleElement[i]=DiagonalMatrix.buildIdentityTimesElementMatrix(i+1, 1/(this.prior.getSD()*this.prior.getSD()));
            }
        innerProductAnswer=new MatrixParameter[LFM.getLoadings().getRowDimension()];
            for (int i = 0; i <innerProductAnswer.length ; i++) {
                innerProductAnswer[i]=new MatrixParameter(null);
                innerProductAnswer[i].setDimensions(i+1,i+1);
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
    }

    private MatrixParameter truncatedInnerProductInPlace(MatrixParameter full, int newRowDimension, MatrixParameter answer){

//        MatrixParameter answer=new MatrixParameter(null);
//        answer.setDimensions(this.getRowDimension(), Right.getRowDimension());
//        System.out.println(answer.getRowDimension());
//        System.out.println(answer.getColumnDimension());

        int p = full.getColumnDimension();
        for (int i = 0; i < newRowDimension; i++) {
            for (int j = 0; j < newRowDimension; j++) {
                double sum = 0;
                for (int k = 0; k < p; k++)
                    sum += full.getParameterValue(i, k) * full.getParameterValue(j,k);
                answer.setParameterValueQuietly(i,j, sum);
            }
        }

        return answer;

    }

    private MatrixParameter truncatedMatrixProductWithTransposedWithVectorInPlace(MatrixParameter Left, double[] Right, int newRowDimension, MatrixParameter answer){

//        MatrixParameter answer=new MatrixParameter(null);
//        answer.setDimensions(this.getRowDimension(), Right.getRowDimension());
//        System.out.println(answer.getRowDimension());
//        System.out.println(answer.getColumnDimension());

        int p = Left.getColumnDimension();
        for (int i = 0; i < newRowDimension; i++) {
            double sum = 0;
                for (int k = 0; k < p; k++)
                    sum += Left.getParameterValue(i, k) * Right[k];
                answer.setParameterValueQuietly(i,0, sum);
            }


        return answer;

    }

    private MatrixParameter getPrecision(int i){
        int size=LFM.getFactorDimension();
        MatrixParameter answer=null;
        if(i<size){
            answer= singleElement[i].add(truncatedInnerProductInPlace(LFM.getFactors(), i + 1, innerProductAnswer[i]).product(LFM.getColumnPrecision().getParameterValue(i, i)));
        }
        else{
             answer= singleElement[size-1].add(innerProductAnswer[size-1].productWithTransposed(innerProductAnswer[size-1]).product(LFM.getColumnPrecision().getParameterValue(i, i)));
            }
        return answer;
    }

    private Matrix getMean(int i, Matrix precision){
//        Matrix factors=null;
        MatrixParameter FxY;
        int size=LFM.getFactorDimension();
        Matrix answer=null;
        double[] scaledDataColumn=LFM.getScaledData().getRowValues(i);
//        Vector dataColumn=null;
//        Vector priorVector=null;
//        Vector temp=null;
//        Matrix data=new Matrix(LFM.getScaledData().getParameterAsMatrix());
        if(i<size){
            FxY=truncatedMatrixProductWithTransposedWithVectorInPlace(LFM.getFactors(), scaledDataColumn, i+1, vectorProductAnswer[i]);
//            dataColumn=new Vector(data.toComponents()[i]);
            try {
                answer=precision.inverse().product(new Matrix(priorMeanVector[i].add(FxY).getParameterAsMatrix()));
            } catch (IllegalDimension illegalDimension) {
                illegalDimension.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            }
        }
        else{
            FxY=truncatedMatrixProductWithTransposedWithVectorInPlace(LFM.getFactors(), scaledDataColumn, size, vectorProductAnswer[size-1]);
//            dataColumn=new Vector(data.toComponents()[i]);
            try {
                answer=precision.inverse().product(new Matrix(priorMeanVector[size-1].add(FxY).getParameterAsMatrix()));
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
            MatrixParameter precision=getPrecision(i);
            Matrix mean=getMean(i, new Matrix(precision.getParameterAsMatrix()));
            draws= MultivariateNormalDistribution.nextMultivariateNormalPrecision(mean.toComponents()[0],precision.getParameterAsMatrix());
            if(i<draws.length){
                while(draws[i]<0){
                    draws= MultivariateNormalDistribution.nextMultivariateNormalPrecision(mean.toComponents()[0],precision.getParameterAsMatrix());
                }
            }
            copy(i, draws);
        }
        LFM.getLoadings().fireParameterChangedEvent();

        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }
}

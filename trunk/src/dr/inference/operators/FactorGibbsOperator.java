package dr.inference.operators;

import dr.evomodel.continuous.LatentFactorModel;
import dr.inference.model.DiagonalMatrix;
import dr.inference.model.MatrixParameter;
import dr.math.matrixAlgebra.IllegalDimension;
import dr.math.matrixAlgebra.Matrix;
import dr.inference.model.Parameter;

/**
 * Created with IntelliJ IDEA.
 * User: max
 * Date: 5/22/14
 * Time: 12:49 PM
 * To change this template use File | Settings | File Templates.
 */
public class FactorGibbsOperator extends SimpleMCMCOperator implements GibbsOperator{
    private static final String FACTOR_GIBBS_OPERATOR="factorGibbsOperator";
    private LatentFactorModel LFM;
    private int numFactors=0;
    private Matrix idMat;
    public FactorGibbsOperator(LatentFactorModel LFM, double weight){
        this.LFM=LFM;
        setWeight(weight);
    }

    private Matrix getPrecision(){
        if(numFactors!=LFM.getFactors().getRowDimension()){
            numFactors=LFM.getFactors().getRowDimension();
            idMat=Matrix.buildIdentity(numFactors);
        }
        Matrix LoadMat=new Matrix(LFM.getLoadings().getParameterAsMatrix());
        Matrix colPrec=new Matrix(LFM.getColumnPrecision().getParameterAsMatrix());
        Matrix answer= null;
        try {
            answer = LoadMat.transpose().product(colPrec).product(LoadMat).add(idMat);
        } catch (IllegalDimension illegalDimension) {
            illegalDimension.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
        return answer;
    }

    private Matrix getMean(){
        Matrix answer=null;
        try {
            answer=getPrecision().inverse().productWithTransposed(new Matrix(LFM.getLoadings().getParameterAsMatrix())).product(new Matrix(LFM.getColumnPrecision().getParameterAsMatrix())).product(LFM.getScaledData());
        } catch (IllegalDimension illegalDimension) {
            illegalDimension.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
        return answer;
    }



    public int getStepCount() {
        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public String getPerformanceSuggestion() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public String getOperatorName() {
        return FACTOR_GIBBS_OPERATOR;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public double doOperation() throws OperatorFailedException {
        Matrix mean=getMean();
        System.out.println(mean.columns());
        System.out.println(mean.rows());
        MatrixParameter meanList=MatrixParameter.parseFromSymmetricDoubleArray(getMean().toComponents());

        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }


}

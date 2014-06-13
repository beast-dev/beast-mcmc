package dr.inference.operators;

import dr.evomodel.continuous.LatentFactorModel;
import dr.inference.distribution.DistributionLikelihood;
import dr.inference.model.DiagonalMatrix;
import dr.inference.model.MatrixParameter;
import dr.math.MathUtils;
import dr.math.distributions.GammaDistribution;

/**
 * Created by max on 6/12/14.
 */
public class LatentFactorModelPrecisionGibbsOperator extends SimpleMCMCOperator implements GibbsOperator {
    private double[] FacXLoad;
    private double[] residual;
    private int factorDimension;
    private LatentFactorModel LFM;
    private GammaDistribution prior;
    private boolean randomScan;

    public LatentFactorModelPrecisionGibbsOperator(LatentFactorModel LFM, DistributionLikelihood prior, double weight, boolean randomScan) {
        setWeight(weight);
        this.LFM=LFM;
        this.prior=(GammaDistribution)prior.getDistribution();
        this.randomScan=randomScan;
        setupParameters();
        FacXLoad=new double[LFM.getFactors().getColumnDimension()];
        residual=new double[LFM.getFactors().getColumnDimension()];
    }

    private void setupParameters(){
        factorDimension=LFM.getFactorDimension();
    }


    private void setPrecision(int i){
        MatrixParameter factors=LFM.getFactors();
        MatrixParameter loadings=LFM.getLoadings();
        DiagonalMatrix precision= (DiagonalMatrix) LFM.getColumnPrecision();
        MatrixParameter data=LFM.getScaledData();
        double shape=prior.getShape()+factors.getColumnDimension()/2;
        double sum;
        for (int j = 0; j <factors.getColumnDimension(); j++) {
            sum=0;
            for (int k = 0; k < factorDimension; k++) {
                sum+=factors.getParameterValue(k,j)*loadings.getParameterValue(k,i);
            }
            FacXLoad[j]=sum;
            residual[j]=data.getParameterValue(i,j)-FacXLoad[j];
        }
        sum=0;
        for (int j = 0; j <factors.getColumnDimension() ; j++) {
            sum+=residual[j]*residual[j];
        }
        double scale=1/(1/prior.getScale()+sum/2);
        double nextPrecision=GammaDistribution.nextGamma(shape, scale);
        precision.getParameter(0).setParameterValueQuietly(i, nextPrecision);
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
        if(LFM.getFactorDimension()!=factorDimension){
            setupParameters();
        }
        if(!randomScan){
        for (int i = 0; i <LFM.getColumnPrecision().getColumnDimension(); i++) {
            setPrecision(i);
        }}
        else{
            int i= MathUtils.nextInt(LFM.getColumnPrecision().getColumnDimension());
            setPrecision(i);
        }
        LFM.getColumnPrecision().getParameter(0).fireParameterChangedEvent();
        
        
        return 0;
    }
}

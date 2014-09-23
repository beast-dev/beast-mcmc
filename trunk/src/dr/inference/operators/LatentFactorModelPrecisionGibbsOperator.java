package dr.inference.operators;

import dr.evomodel.continuous.LatentFactorModel;
import dr.inference.distribution.DistributionLikelihood;
import dr.inference.model.DiagonalMatrix;
import dr.inference.model.MaskedParameter;
import dr.inference.model.MatrixParameter;
import dr.math.MathUtils;
import dr.math.distributions.GammaDistribution;

/**
 * Created by max on 6/12/14.
 */
public class LatentFactorModelPrecisionGibbsOperator extends SimpleMCMCOperator implements GibbsOperator {
//    private double[] FacXLoad;
//    private double[] residual;
    private LatentFactorModel LFM;
    private GammaDistribution prior;
    private boolean randomScan;
    private double shape;

    public LatentFactorModelPrecisionGibbsOperator(LatentFactorModel LFM, DistributionLikelihood prior, double weight, boolean randomScan) {
        setWeight(weight);
        this.LFM=LFM;
        this.prior=(GammaDistribution)prior.getDistribution();
        this.randomScan=randomScan;

//        FacXLoad=new double[LFM.getFactors().getColumnDimension()];
//        residual=new double[LFM.getFactors().getColumnDimension()];
        shape=this.prior.getShape()+LFM.getFactors().getColumnDimension()*.5;
    }




    private void setPrecision(int i){
        MatrixParameter factors=LFM.getFactors();
        MatrixParameter loadings=LFM.getLoadings();
        DiagonalMatrix precision= (DiagonalMatrix) LFM.getColumnPrecision();
        MatrixParameter data=LFM.getScaledData();
        double di=0;
        for (int j = 0; j <factors.getColumnDimension(); j++) {
            double sum=0;
            for (int k = 0; k < factors.getRowDimension(); k++) {
                sum+=factors.getParameterValue(k,j)*loadings.getParameterValue(k,i);
            }
            double temp=data.getParameterValue(i,j)-sum;
            di+=temp*temp;
//            FacXLoad[j]=sum;
//            residual[j]=data.getParameterValue(i,j)-FacXLoad[j];
//        }
//        double sum=0;
//        for (int j = 0; j <factors.getColumnDimension() ; j++) {
//            sum+=residual[j]*residual[j];
        }
        double scale=1.0/(1.0/prior.getScale()+di*.5);
        double nextPrecision=GammaDistribution.nextGamma(shape, scale);
        precision.setParameterValueQuietly(i, nextPrecision);
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

        if(!randomScan) for (int i = 0; i < LFM.getColumnPrecision().getColumnDimension(); i++) {
            if (!(LFM.getColumnPrecision().getParameter(0) instanceof MaskedParameter))
             setPrecision(i);
            else{
                if((((MaskedParameter) LFM.getColumnPrecision().getParameter(0)).getParameterMaskValue(i))!=0)
                    setPrecision(i);
            }
        }
        else{
            int i= MathUtils.nextInt(LFM.getColumnPrecision().getColumnDimension());
            if ((LFM.getColumnPrecision().getParameter(0) instanceof MaskedParameter)){
            while((((MaskedParameter) LFM.getColumnPrecision().getParameter(0)).getParameterMaskValue(i))==0)
                i= MathUtils.nextInt(LFM.getColumnPrecision().getColumnDimension());     }
            setPrecision(i);
        }
        LFM.getColumnPrecision().getParameter(0).fireParameterChangedEvent();
        
        
        return 0;
    }
}

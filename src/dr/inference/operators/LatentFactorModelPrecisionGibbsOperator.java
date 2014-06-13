package dr.inference.operators;

import dr.evomodel.continuous.LatentFactorModel;
import dr.inference.distribution.DistributionLikelihood;
import dr.math.distributions.GammaDistribution;

/**
 * Created by max on 6/12/14.
 */
public class LatentFactorModelPrecisionGibbsOperator extends SimpleMCMCOperator implements GibbsOperator {
    private double[] FacXLoad;
    private double[] residual;
    private LatentFactorModel LFM;
    private GammaDistribution prior;
    private boolean randomScan;

    public LatentFactorModelPrecisionGibbsOperator(LatentFactorModel LFM, DistributionLikelihood prior, double weight, boolean randomScan) {
        setWeight(weight);
        this.LFM=LFM;
        this.prior=(GammaDistribution)prior.getDistribution();
        this.randomScan=randomScan;
        setupParameters();
    }

    private void setupParameters(){
        FacXLoad=new double[LFM.getFactorDimension()];
        residual=new double[LFM.getFactorDimension()];
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
        return 0;
    }
}

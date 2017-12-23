package dr.inference.operators;

import dr.inference.distribution.DistributionLikelihood;
import dr.inference.distribution.IndependentInverseGammaDistributionModel;
import dr.inference.model.Parameter;
import dr.math.distributions.GammaDistribution;
import dr.math.distributions.InverseGammaDistribution;

public class ShrinkageAugmentedGibbsOperator extends SimpleMCMCOperator implements GibbsOperator {

    Parameter localShape;
    Parameter globalShape;

    Parameter localAugmented;
    Parameter globalAugmented;
    Parameter local;
    Parameter global;

    double localAugmentedShape;
    double localAugmentedRate;
    double globalAugmentedShape;
    double globalAugmentedRate;

    double pathParameter;


    public ShrinkageAugmentedGibbsOperator(double weight, DistributionLikelihood localAugmentedPrior, DistributionLikelihood globalAugmentedPrior,
                                           IndependentInverseGammaDistributionModel localPrior,
                                           IndependentInverseGammaDistributionModel globalPrior){
        setWeight(weight);
        GammaDistribution localAugmentedPriorGamma = (GammaDistribution) localAugmentedPrior.getDistribution();
        GammaDistribution globalAugmentedPriorGamma = (GammaDistribution) globalAugmentedPrior.getDistribution();

        this.localShape = localPrior.getShape();
        this.globalShape = globalPrior.getShape();

        this.local = localPrior.getData();
        this.global = globalPrior.getData();
        this.localAugmented = localPrior.getScale();
        this.globalAugmented = globalPrior.getScale();

        localAugmentedShape = localAugmentedPriorGamma.getShape();
        localAugmentedRate = 1 / localAugmentedPriorGamma.getScale();
        globalAugmentedShape = globalAugmentedPriorGamma.getShape();
        globalAugmentedRate = 1 / globalAugmentedPriorGamma.getScale();
        this.pathParameter = 1;
    }

    public int getStepCount() {
        return 0;
    }

    @Override
    public String getPerformanceSuggestion() {
        return null;
    }

    @Override
    public String getOperatorName() {
        return "ShrinkageAugmentedGibbsOperator";
    }

    @Override
    public double doOperation() {
        for (int i = 0; i < local.getDimension(); i++) {
            double shape = localAugmentedShape +  localShape.getParameterValue(i);
            double rate = localAugmentedRate +  1 / local.getParameterValue(i);
//            double scale = 1 / localAugmentedRate + pathParameter * local.getParameterValue(i);
            localAugmented.setParameterValueQuietly(i, GammaDistribution.nextGamma(shape, 1 / rate));
        }
        localAugmented.fireParameterChangedEvent();

        double shape = globalAugmentedShape +  globalShape.getParameterValue(0);
        double rate = globalAugmentedRate +  1 / global.getParameterValue(0);
//        double scale = 1 / globalAugmentedRate + pathParameter *  global.getParameterValue(0);
        globalAugmented.setParameterValue(0, GammaDistribution.nextGamma(shape, 1 / rate));

        return 0;
    }

}

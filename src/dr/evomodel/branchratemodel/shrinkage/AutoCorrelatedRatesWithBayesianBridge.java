package dr.evomodel.branchratemodel.shrinkage;

import dr.evomodel.branchratemodel.AutoCorrelatedBranchRatesDistribution;
import dr.inference.distribution.shrinkage.BayesianBridgeDistributionModel;
import dr.inference.distribution.shrinkage.BayesianBridgeStatisticsProvider;
import dr.inference.model.Parameter;

/**
 * @author Marc A. Suchard
 */

@Deprecated
public class AutoCorrelatedRatesWithBayesianBridge implements BayesianBridgeStatisticsProvider {

    private final AutoCorrelatedBranchRatesDistribution rates;
    private final BayesianBridgeDistributionModel prior;

    public AutoCorrelatedRatesWithBayesianBridge(AutoCorrelatedBranchRatesDistribution rates,
                                                 BayesianBridgeDistributionModel prior) {
        this.rates = rates;
        this.prior = prior;
    }

    @Override
    public double getCoefficient(int i) {
        return rates.getIncrement(i);
    }

    @Override
    public Parameter getGlobalScale() {
        return prior.getGlobalScale();
    }

    @Override
    public Parameter getLocalScale() {
        return prior.getLocalScale();
    }

    @Override
    public Parameter getExponent() {
        return prior.getExponent();
    }

    @Override
    public Parameter getSlabWidth() {
        return prior.getSlabWidth();
    }

    @Override
    public int getDimension() {
        return rates.getDimension();
    }
}

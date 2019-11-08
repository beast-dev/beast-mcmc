package dr.inference.distribution.shrinkage;

import dr.inference.model.Parameter;

/**
 * @author Marc A. Suchard
 */
public interface BayesianBridgeStatisticsProvider {

    double getCoefficient(int i);

    Parameter getGlobalScale();

    Parameter getLocalScale();

    Parameter getExponent();

    int getDimension();
}

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

    Parameter getSlabWidth();

    int getDimension();

    static boolean equivalent(BayesianBridgeStatisticsProvider lhs,
                              BayesianBridgeStatisticsProvider rhs) {
        return
                lhs.getDimension() == rhs.getDimension()  &&
                lhs.getExponent() == rhs.getExponent() &&
                lhs.getGlobalScale() == rhs.getGlobalScale() &&
                lhs.getLocalScale() == rhs.getLocalScale() &&
                lhs.getSlabWidth() == rhs.getSlabWidth();
    }
}

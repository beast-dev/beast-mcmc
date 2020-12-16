package dr.evomodel.treelikelihood.thorneytreelikelihood;

import dr.inference.model.*;
import org.apache.commons.math.special.Gamma;
import org.apache.commons.math.util.FastMath;
import org.apache.commons.math.util.MathUtils;

public class StrictClockBranchLengthLikelihoodDelegate extends AbstractModel implements ThorneyBranchLengthLikelihoodDelegate {
    private Parameter mutationRate;
    private final double scale;

    public StrictClockBranchLengthLikelihoodDelegate(String name, Parameter mutationRate, double scale){
        super(name);
        this.mutationRate = mutationRate;
        this.scale = scale;
        addVariable(mutationRate);

    }

    @Override
    public double getLogLikelihood(double mutations, double time) {
        return SaddlePointExpansion.logPoissonProbability(time*mutationRate.getValue(0)*scale, (int) Math.round(mutations));
    }


    @Override
    public double getGradientWrtTime(double mutations, double time) { // TODO: better chain rule handling
        return SaddlePointExpansion.logPoissonMeanDerivative(time * mutationRate.getValue(0) * scale, (int) Math.round(mutations)) * mutationRate.getValue(0) * scale;
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {

    }

    /**
     * This method is called whenever a parameter is changed.
     * <p/>
     * It is strongly recommended that the model component sets a "dirty" flag and does no
     * further calculations. Recalculation is typically done when the model component is asked for
     * some information that requires them. This mechanism is 'lazy' so that this method
     * can be safely called multiple times with minimal computational cost.
     *
     * @param variable
     * @param index
     * @param type
     */
    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        fireModelChanged();
    }

    /**
     * Additional state information, outside of the sub-model is stored by this call.
     */
    @Override
    protected void storeState() {

    }

    /**
     * After this call the model is guaranteed to have returned its extra state information to
     * the values coinciding with the last storeState call.
     * Sub-models are handled automatically and do not need to be considered in this method.
     */
    @Override
    protected void restoreState() {

    }

    /**
     * This call specifies that the current state is accept. Most models will not need to do anything.
     * Sub-models are handled automatically and do not need to be considered in this method.
     */
    @Override
    protected void acceptState() {

    }


}

// Grabbed some stuff from Commons Maths as it is not public
// This code is under the Apache License 2.0

final class SaddlePointExpansion {
    private static final double HALF_LOG_2_PI = 0.5D * FastMath.log(6.283185307179586D);
    private static final double[] EXACT_STIRLING_ERRORS = new double[]{0.0D, 0.15342640972002736D, 0.08106146679532726D, 0.05481412105191765D, 0.0413406959554093D, 0.03316287351993629D, 0.02767792568499834D, 0.023746163656297496D, 0.020790672103765093D, 0.018488450532673187D, 0.016644691189821193D, 0.015134973221917378D, 0.013876128823070748D, 0.012810465242920227D, 0.01189670994589177D, 0.011104559758206917D, 0.010411265261972096D, 0.009799416126158804D, 0.009255462182712733D, 0.008768700134139386D, 0.00833056343336287D, 0.00793411456431402D, 0.007573675487951841D, 0.007244554301320383D, 0.00694284010720953D, 0.006665247032707682D, 0.006408994188004207D, 0.006171712263039458D, 0.0059513701127588475D, 0.0057462165130101155D, 0.005554733551962801D};

    private SaddlePointExpansion() {
    }

    static double getStirlingError(double z) {
        double ret;
        double z2;
        if (z < 15.0D) {
            z2 = 2.0D * z;
            if (FastMath.floor(z2) == z2) {
                ret = EXACT_STIRLING_ERRORS[(int)z2];
            } else {
                ret = Gamma.logGamma(z + 1.0D) - (z + 0.5D) * FastMath.log(z) + z - HALF_LOG_2_PI;
            }
        } else {
            z2 = z * z;
            ret = (0.08333333333333333D - (0.002777777777777778D - (7.936507936507937E-4D - (5.952380952380953E-4D - 8.417508417508417E-4D / z2) / z2) / z2) / z2) / z;
        }

        return ret;
    }

    static double getDeviancePart(double x, double mu) {
        double ret;
        if (FastMath.abs(x - mu) < 0.1D * (x + mu)) {
            double d = x - mu;
            double v = d / (x + mu);
            double s1 = v * d;
            double s = 0.0D / 0.0;
            double ej = 2.0D * x * v;
            v *= v;

            for(int j = 1; s1 != s; ++j) {
                s = s1;
                ej *= v;
                s1 += ej / (double)(j * 2 + 1);
            }

            ret = s1;
        } else {
            ret = x * FastMath.log(x / mu) + mu - x;
        }

        return ret;
    }

    static double logBinomialProbability(int x, int n, double p, double q) {
        double ret;
        if (x == 0) {
            if (p < 0.1D) {
                ret = -getDeviancePart((double)n, (double)n * q) - (double)n * p;
            } else {
                ret = (double)n * FastMath.log(q);
            }
        } else if (x == n) {
            if (q < 0.1D) {
                ret = -getDeviancePart((double)n, (double)n * p) - (double)n * q;
            } else {
                ret = (double)n * FastMath.log(p);
            }
        } else {
            ret = getStirlingError((double)n) - getStirlingError((double)x) - getStirlingError((double)(n - x)) - getDeviancePart((double)x, (double)n * p) - getDeviancePart((double)(n - x), (double)n * q);
            double f = 6.283185307179586D * (double)x * (double)(n - x) / (double)n;
            ret += -0.5D * FastMath.log(f);
        }

        return ret;
    }

    static public double logPoissonProbability(double mean,int x) {
//        double ret;
//        if (x >= 0 && x != Integer.MAX_VALUE) {
//            if (x == 0) {
//                ret = FastMath.exp(-mean);
//            } else {
//                ret = FastMath.exp(-getStirlingError((double)x) - getDeviancePart((double)x, mean)) / FastMath.sqrt(6.283185307179586D * (double)x);
//            }
//        } else {
//            ret = 0.0D;
//        }
//        ret = Math.log(ret);

        double result;
        if (x >= 0 && x != Integer.MAX_VALUE) {
            if (x == 0) {
                result = -mean;
            } else {
                result = -getStirlingError((double)x) - getDeviancePart((double)x, mean) - Math.log(MathUtils.TWO_PI * (double)x) * 0.5;
            }
        } else {
            result = Double.NEGATIVE_INFINITY;
        }

        return result;
    }

    static public double logPoissonMeanDerivative(double mean, int x) {
        final double result = x == 0 ? -1.0 : FastMath.log(((double) x) / mean) - 1.0;

        return result;
    }
}

package dr.evomodel.operators;

import dr.inference.hmc.PrecisionColumnProvider;
import dr.inference.operators.hmc.MinimumTravelInformation;

/**
 * @author Marc A. Suchard
 */
public class NativeZigZagWrapper {

    private final int instanceNumber;

    public NativeZigZagWrapper(int dimension,
                               PrecisionColumnProvider columnProvider,
                               double[] mask,
                               double[] observed) {
        this.instanceNumber = NativeZigZag.INSTANCE.createInstance(dimension, columnProvider, mask, observed);
    }

    public void operate(double[] position,
                        double[] velocity,
                        double[] action,
                        double[] gradient,
                        double[] moment,
                        double time) {
        NativeZigZag.INSTANCE.operate(instanceNumber, position, velocity, action, gradient, moment, time);
    }

    public MinimumTravelInformation testGetNextEvent(double[] position,
                                                     double[] velocity,
                                                     double[] action,
                                                     double[] gradient,
                                                     double[] moment) {
        return NativeZigZag.INSTANCE.getNextEvent(instanceNumber, position, velocity, action, gradient, moment);
    }
}

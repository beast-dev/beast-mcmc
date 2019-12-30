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

    public void operate(PrecisionColumnProvider columnProvider,
                        double[] position,
                        double[] velocity,
                        double[] action,
                        double[] gradient,
                        double[] moment,
                        double time) {
        NativeZigZag.INSTANCE.operate(instanceNumber, columnProvider, position, velocity, action, gradient, moment, time);
    }

    public MinimumTravelInformation getNextEvent(double[] position,
                                                 double[] velocity,
                                                 double[] action,
                                                 double[] gradient,
                                                 double[] moment) {
        return NativeZigZag.INSTANCE.getNextEvent(instanceNumber, position, velocity, action, gradient, moment);
    }

    public int enterCriticalRegion(
            double[] position,
            double[] velocity,
            double[] action,
            double[] gradient,
            double[] momentum) {
        return NativeZigZag.INSTANCE.enterCriticalRegion(instanceNumber, position, velocity, action, gradient, momentum);
    }

    public int exitCriticalRegion() {
        return NativeZigZag.INSTANCE.exitCriticalRegion(instanceNumber);
    }

    public boolean inCriticalRegion() {
        return NativeZigZag.INSTANCE.inCriticalRegion(instanceNumber);
    }

    public MinimumTravelInformation getNextEventInCriticalRegion() {
        return NativeZigZag.INSTANCE.getNextEventInCriticalRegion(instanceNumber);
    }

    public void innerBounce(double[] position,
                            double[] velocity,
                            double[] action,
                            double[] gradient,
                            double[] momentum,
                            double eventTime, int eventIndex, int eventType) {
        NativeZigZag.INSTANCE.innerBounce(instanceNumber, position, velocity, action, gradient, momentum,
                eventTime, eventIndex, eventType);
    }

    public void innerBounceCriticalRegion(double eventTime, int eventIndex, int eventType) {
        NativeZigZag.INSTANCE.innerBounceCriticalRegion(instanceNumber, eventTime, eventIndex, eventType);
    }

    public void updateDynamics(double[] position,
                               double[] velocity,
                               double[] action,
                               double[] gradient,
                               double[] momentum,
                               double[] column,
                               double eventTime, int eventIndex, int eventType) {
        NativeZigZag.INSTANCE.updateDynamics(instanceNumber, position, velocity, action, gradient, momentum,
                column, eventTime, eventIndex, eventType);
    }
}

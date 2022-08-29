package dr.evomodel.operators;

import dr.inference.hmc.PrecisionColumnProvider;
import dr.inference.operators.hmc.MinimumTravelInformation;
import dr.inference.operators.hmc.MinimumTravelInformationBinary;

/**
 * @author Marc A. Suchard
 */
public class NativeZigZagWrapper {

    private final int instanceNumber;

    public NativeZigZagWrapper(int dimension,
                               NativeZigZagOptions options,
                               double[] mask,
                               double[] observed,
                               double[] parameterSign) {
        this.instanceNumber = NativeZigZag.INSTANCE.createInstance(dimension, options, mask, observed, parameterSign);
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

    public MinimumTravelInformationBinary getNextReversibleEvent(double[] position,
                                                                 double[] velocity,
                                                                 double[] action,
                                                                 double[] gradient,
                                                                 double[] momentum) {
        return NativeZigZag.INSTANCE.getNextEvent(instanceNumber, position, velocity, action, gradient, momentum);
    }

    public MinimumTravelInformation getNextIrreversibleEvent(double[] position,
                                                             double[] velocity,
                                                             double[] action,
                                                             double[] gradient) {
        throw new RuntimeException("not implemented yet");
    }

    public void updateReversibleDynamics(double[] position,
                                         double[] velocity,
                                         double[] action,
                                         double[] gradient,
                                         double[] momentum,
                                         double[] column,
                                         double eventTime, int eventIndex, int eventType) {
        NativeZigZag.INSTANCE.updateDynamics(instanceNumber, position, velocity, action, gradient, momentum,
                column, eventTime, eventIndex, eventType);
    }

    public void updateIrreversibleDynamics(double[] position,
                                           double[] velocity,
                                           double[] action,
                                           double[] gradient,
                                           double[] column,
                                           double eventTime, int eventIndex, int eventType) {
        NativeZigZag.INSTANCE.updateDynamics(instanceNumber, position, velocity, action, gradient, null,
                column, eventTime, eventIndex, eventType);
    }

    @SuppressWarnings("unused")
    public int enterCriticalRegion(
            double[] position,
            double[] velocity,
            double[] action,
            double[] gradient,
            double[] momentum) {
        return NativeZigZag.INSTANCE.enterCriticalRegion(instanceNumber, position, velocity, action, gradient, momentum);
    }

    @SuppressWarnings("unused")
    public int exitCriticalRegion() {
        return NativeZigZag.INSTANCE.exitCriticalRegion(instanceNumber);
    }

    @SuppressWarnings("unused")
    public boolean inCriticalRegion() {
        return NativeZigZag.INSTANCE.inCriticalRegion(instanceNumber);
    }

    @SuppressWarnings("unused")
    public MinimumTravelInformation getNextEventInCriticalRegion() {
        return NativeZigZag.INSTANCE.getNextEventInCriticalRegion(instanceNumber);
    }

    @SuppressWarnings("unused")
    public void innerBounce(double[] position,
                            double[] velocity,
                            double[] action,
                            double[] gradient,
                            double[] momentum,
                            double eventTime, int eventIndex, int eventType) {
        NativeZigZag.INSTANCE.innerBounce(instanceNumber, position, velocity, action, gradient, momentum,
                eventTime, eventIndex, eventType);
    }

    @SuppressWarnings("unused")
    public void innerBounceCriticalRegion(double eventTime, int eventIndex, int eventType) {
        NativeZigZag.INSTANCE.innerBounceCriticalRegion(instanceNumber, eventTime, eventIndex, eventType);
    }
}

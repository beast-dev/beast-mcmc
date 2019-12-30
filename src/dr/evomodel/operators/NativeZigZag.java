package dr.evomodel.operators;

import dr.inference.hmc.PrecisionColumnProvider;
import dr.inference.operators.hmc.MinimumTravelInformation;

class NativeZigZag {

    static NativeZigZag INSTANCE;
    private static int instanceNumber = 0;

    private NativeZigZag() { }

    int createInstance(int dimension,
                       PrecisionColumnProvider columnProvider,
                       double[] mask,
                       double[] observed) {

        double[] column = columnProvider.getColumn(0);
        if (dimension != column.length || (mask != null && dimension != mask.length) ||
                (observed != null && dimension != observed.length)) {
            throw new RuntimeException("Invalid dimensions");
        }

        int result = create(dimension, columnProvider, mask, observed);
        if (result < 0) {
            throw new RuntimeException("Unable to create instance");
        }

        return instanceNumber++;
    }

    private native int create(int dimension,
                              PrecisionColumnProvider columnProvider,
                              double[] mask,
                              double[] observed);

    native int operate(int instanceNumber,
                       PrecisionColumnProvider columnProvider,
                       double[] position,
                       double[] velocity,
                       double[] action,
                       double[] gradient,
                       double[] momentum,
                       double time);

    native MinimumTravelInformation getNextEvent(int instanceNumber,
                                                 double[] position,
                                                 double[] velocity,
                                                 double[] action,
                                                 double[] gradient,
                                                 double[] momentum);

    native int enterCriticalRegion(int instanceNumber,
                                   double[] position,
                                   double[] velocity,
                                   double[] action,
                                   double[] gradient,
                                   double[] momentum);

    native int exitCriticalRegion(int instanceNumber);

    native boolean inCriticalRegion(int instanceNumber);

    native MinimumTravelInformation getNextEventInCriticalRegion(int instanceNumber);

    native void innerBounce(int instanceNumber,
                            double[] position,
                            double[] velocity,
                            double[] action,
                            double[] gradient,
                            double[] momentum,
                            double eventTime, int eventIndex, int eventType);

    native void innerBounceCriticalRegion(int instanceNumber,
                                          double eventTime, int eventIndex, int eventType);

    native void updateDynamics(int instanceNumber,
                               double[] position,
                               double[] velocity,
                               double[] action,
                               double[] gradient,
                               double[] momentum,
                               double[] column,
                               double eventTime, int eventIndex, int eventType);

    static {
        System.loadLibrary("zig_zag");
        INSTANCE = new NativeZigZag();
    }


}

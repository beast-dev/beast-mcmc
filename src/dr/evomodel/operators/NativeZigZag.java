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
                                                 double[] momentum,
                                                 int type,
                                                 int index,
                                                 double time);

    static {
        System.loadLibrary("zig_zag");
        INSTANCE = new NativeZigZag();
    }

}

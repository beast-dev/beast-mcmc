package dr.evomodel.operators;

import dr.inference.hmc.PrecisionColumnProvider;
import dr.inference.operators.hmc.MinimumTravelInformation;

public class NativeZigZag {

    public enum Flag {
        PRECISION_DOUBLE(1L << 1, "double precision computation"),
        PRECISION_SINGLE(1L << 2, "single precision computation"),
        FRAMEWORK_TBB(1L << 3, "TBB multi-core threading"),
        FRAMEWORK_OPENCL(1L << 4, "use OpenCL implementation with CPU or GPU resources"),
        SIMD_SSE(1L << 7, "use SSE SIMD vectorization"),
        SIMD_AVX(1L << 8, "use AVX SIMD vectorization"),
        SIMD_AVX512(1L << 9, "use AVX-512 SIMD vectorization");

        private final long mask;
        private final String meaning;

        Flag(long mask, String meaning) {
            this.mask = mask;
            this.meaning = meaning;
        }

        @SuppressWarnings("unused")
        public String getMeaning() { return meaning; }

        public long getMask() {
            return this.mask;
        }

        public boolean isSet(long flag) {
            return (flag & this.mask) != 0L;
        }

        public static String toString(long flag) {
            StringBuilder sb = new StringBuilder();
            for (Flag mask : values()) {
                if (mask.isSet(flag)) {
                    sb.append(" ").append(mask.name());
                }
            }

            return sb.toString();
        }
    }

    static NativeZigZag INSTANCE;
    private static int instanceNumber = 0;

    private NativeZigZag() { }

    int createInstance(int dimension,
                       NativeZigZagOptions options,
                       double[] mask,
                       double[] observed) {


        if ((mask != null && dimension != mask.length) ||
                (observed != null && dimension != observed.length)) {
            throw new RuntimeException("Invalid dimensions");
        }

        int result = create(dimension, options, mask, observed);
        if (result < 0) {
            throw new RuntimeException("Unable to create instance");
        }

        return instanceNumber++;
    }

    private native int create(int dimension,
                              NativeZigZagOptions options,
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

    native MinimumTravelInformation getNextEventIrreversible(int instanceNumber,
                                                             double[] position,
                                                             double[] velocity,
                                                             double[] action,
                                                             double[] gradient);
    static {
        System.loadLibrary("zig_zag");
        INSTANCE = new NativeZigZag();
    }
}

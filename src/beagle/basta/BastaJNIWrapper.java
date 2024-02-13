package beagle.basta;

public class BastaJNIWrapper {

    private static final String LIBRARY_NAME = getPlatformSpecificLibraryName();

    private BastaJNIWrapper() { }

    public native int allocateCoalescentBuffers(int instance,
                                                int bufferCount,
                                                int maxCoalescentIntervalCount); // TODO buffers have different sizes

    public native int getBastaBuffer(int instance,
                                     int index,
                                     double[] buffer);

    public native int updateBastaPartials(int instance,
                                          final int[] operations,
                                          int operationCount,
                                          final int[] intervals,
                                          int intervalCount,
                                          int populationSizeIndex,
                                          int coalescentProbabilityIndex);

    public native int accumulateBastaPartials(int instance,
                                              final int[] operations,
                                              int operationCount,
                                              final int[] intervals,
                                              int intervalCount,
                                              final double[] intervalLengths,
                                              int populationSizesIndex,
                                              int coalescentProbabilityIndex,
                                              double[] result);

    private static String getPlatformSpecificLibraryName() {
        String osName = System.getProperty("os.name").toLowerCase();
        String osArch = System.getProperty("os.arch").toLowerCase();
        if (osName.startsWith("windows")) {
            if (osArch.equals("x86") || osArch.equals("i386")) return "hmsbeagle-basta32";
            if (osArch.startsWith("amd64") || osArch.startsWith("x86_64")) return "hmsbeagle-basta64";
        }
        return "hmsbeagle-jni-basta";
    }

    public static void loadBastaLibrary() throws UnsatisfiedLinkError {
        String path = "";
        if (System.getProperty("beagle.library.path") != null) {
            path = System.getProperty("beagle.library.path");
            if (path.length() > 0 && !path.endsWith("/")) {
                path += "/";
            }
        }

        System.loadLibrary(path + LIBRARY_NAME);
        INSTANCE = new BastaJNIWrapper();
    }

    public static BastaJNIWrapper INSTANCE;
}

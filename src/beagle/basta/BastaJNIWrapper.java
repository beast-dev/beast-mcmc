package beagle.basta;

public class BastaJNIWrapper {

    private static final String LIBRARY_NAME = getPlatformSpecificLibraryName();

    private BastaJNIWrapper() { }

    public native int updateBastaPartials(int instance,
                                          final int[] operations,
                                          int operationCount,
                                          int populationSizeIndex);

    public native int accumulateBastaPartials(int instance,
                                              final int[] operations,
                                              int operationCount,
                                              final int[] segments,
                                              int segmentCount);

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

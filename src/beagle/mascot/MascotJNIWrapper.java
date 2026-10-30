/*
 * MascotJNIWrapper.java
 *
 * Copyright (c) 2002-2026 Alexei Drummond, Andrew Rambaut and Marc Suchard
 */

package beagle.mascot;

import java.io.File;

/**
 * The MASCOT component of beagle-lib (libhmsbeagle/mascot/mascot.h, API v1) through its JNI library
 * {@code hmsbeagle-jni-mascot}. The library is loaded from the directory named by the {@code mascot.library.path}
 * system property if it is set, otherwise from {@code java.library.path}. Constants mirror mascot.h.
 *
 * @author Marc A. Suchard
 */
public final class MascotJNIWrapper {

    public static final String LIBRARY_NAME = "hmsbeagle-jni-mascot";
    public static final String LIBRARY_PATH_PROPERTY = "mascot.library.path";
    public static final int API_VERSION_MAJOR = 1;

    public static final long FP_STRICT = 1L;
    public static final long FP_FMA = 1L << 1;
    public static final long FORCE_KERNEL_A = 1L << 2;
    public static final long FORCE_KERNEL_B = 1L << 3;

    public static final int OUT_LOGL = 1;
    public static final int OUT_THETA_GRADIENT = 1 << 1;
    public static final int OUT_CLOCK_GRADIENT = 1 << 2;
    public static final int OUT_ANCESTRAL_SCORES = 1 << 3;

    public static final int OK = 0;
    public static final int ERROR_ARGUMENT = -1;
    public static final int ERROR_COALESCENT_RATE = -2;
    public static final int ERROR_PROBABILITIES = -3;
    public static final int ERROR_UNSUPPORTED = -4;
    public static final int ERROR_OUT_OF_MEMORY = -5;

    public static final int STEPS_FIXED = 0;
    public static final int STEPS_ADAPTIVE = 1;

    private static final String LOAD_FAILURE = load();

    private MascotJNIWrapper() {
    }

    private static String load() {
        try {
            String directory = System.getProperty(LIBRARY_PATH_PROPERTY);
            if (directory != null) {
                System.load(new File(directory, System.mapLibraryName(LIBRARY_NAME)).getAbsolutePath());
            } else {
                System.loadLibrary(LIBRARY_NAME);
            }
            return null;
        } catch (UnsatisfiedLinkError | SecurityException e) {
            return e.getMessage();
        }
    }

    /**
     * @return null if the library is loaded and implements API version {@value #API_VERSION_MAJOR}.x, otherwise why
     * it cannot be used
     */
    public static String getUnavailableReason() {
        if (LOAD_FAILURE != null) {
            return "the " + LIBRARY_NAME + " library could not be loaded (" + LOAD_FAILURE + ")";
        }
        int[] version = new int[2];
        if (getApiVersion(version) != OK || version[0] != API_VERSION_MAJOR) {
            return "the " + LIBRARY_NAME + " library implements API " + version[0] + "." + version[1] +
                    ", not " + API_VERSION_MAJOR + ".x";
        }
        return null;
    }

    public static native int getApiVersion(int[] majorMinor);

    /** @return an instance id >= 0, or a negative status */
    public static native int createInstance(int stateCount, int lineageCount, int maxEpochCount, long flags);

    public static native int finalizeInstance(int instance);

    /** @return "kernel simdWidth fpMode isa", e.g. "A 0 fma neon", or null for an unknown instance */
    public static native String getInstanceDetails(int instance);

    /** @param tipPartials lineageCount x stateCount; rows of internal nodes are ignored */
    public static native int setTipPartials(int instance, double[] tipPartials);

    /**
     * One evaluation (mascotEvaluate): forward only when {@code outputs == OUT_LOGL}, otherwise forward, tape and
     * adjoint. Arrays are copied in and out; {@code branchRates} and unrequested outputs may be null.
     */
    public static native int evaluate(int instance,
                                      int intervalCount, double initialTime, int initialSampleNode,
                                      double[] intervalLengths, int[] eventTypes, int[] sampleNodes,
                                      int[] coalescentChild1, int[] coalescentChild2, int[] coalescentParent,
                                      int epochCount, double[] boundaries, double[] migrationMatrices,
                                      double[] inversePopulations, double[] branchRates,
                                      int stepRuleKind, double maxStep, double tau, double hMin,
                                      int maxStepsPerInterval,
                                      boolean checkProbabilities, int checkpointReadSet, int checkpointWriteSet,
                                      int outputs,
                                      double[] logLikelihood, double[] thetaGradient, double[] clockGradient,
                                      double[] ancestralStateScores);
}

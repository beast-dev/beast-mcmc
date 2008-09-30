package dr.evomodel.treelikelihood;

import java.util.logging.Logger;

/**
 * @author Marc Suchard
 * @author Dat Hunyh
 *
 */

public class GPUMemoryLikelihoodCore extends NativeMemoryLikelihoodCore {

    /**
     * Constructor
     *
     * @param stateCount number of states
     */
    public GPUMemoryLikelihoodCore(int stateCount) {
        super(stateCount);
        StringBuffer sb = new StringBuffer();
        sb.append("Constructing GPU likelihood core:\n");
        sb.append("\tGPU Name: ");
        sb.append("\tOther info here");
        sb.append("If you publish results using this core, please reference Hunyh and Suchard (in preparation)");
        Logger.getLogger("dr.evomodel.treelikelihood").info(sb.toString());
    }

    /* GPU memory handing functions */

    private native long allocateNativeMemoryArray(int length);

    private native long allocateNativeIntMemoryArray(int length);

    private native void freeNativeMemoryArray(long nativePtr);

    private native void setNativeMemoryArray(double[] fromArray, int fromOffset,
                                             long toNativePtr, int toOffset, int length);

    private native void setNativeMemoryArray(int[] fromArray, int fromOffset,
                                             long toNativePtr, int toOffset, int length);

    private native void getNativeMemoryArray(long fromNativePtr, int fromOffset,
                                             double[] toArray, int toOffset, int length);

    private native void getNativeMemoryArray(long fromNativePtr, int fromOffset,
                                             int[] toArray, int toOffset, int length);


    /* GPU peeling functions */

    public native void nativeIntegratePartials(long ptrPartials, double[] proportions,
                                               int patternCount, int matrixCount,
                                               double[] outPartials,
                                               int stateCount);

    protected native void nativePartialsPartialsPruning(long ptrPartials1, long ptrMatrices1,
                                                        long ptrPartials2, long ptrMatrices2,
                                                        int patternCount, int matrixCount,
                                                        long ptrPartials3,
                                                        int stateCount);

    protected native void nativeStatesPartialsPruning(long ptrStates1, long ptrMatrices1,
                                                      long ptrPartials2, long ptrMatrices2,
                                                      int patternCount, int matrixCount,
                                                      long ptrPartials3,
                                                      int stateCount);

    protected native void nativeStatesStatesPruning(long ptrStates1, long ptrMatrices1,
                                                    long ptrStates2, long ptrMatrices2,
                                                    int patternCount, int matrixCount,
                                                    long ptrPartials3,
                                                    int stateCount);

    /* Library loading routines */

    public static boolean isAvailable() {
        return isNativeAvailable;
    }

    public static boolean isGPUCompatiable() {
        // todo Test for NVIDIA, etc.
        return true;
    }

    private static boolean isNativeAvailable = false;

    static {
        try {
            System.loadLibrary("NativeMemoryLikelihoodCore");
            isNativeAvailable = true;
        } catch (UnsatisfiedLinkError e) {
        }
    }
}

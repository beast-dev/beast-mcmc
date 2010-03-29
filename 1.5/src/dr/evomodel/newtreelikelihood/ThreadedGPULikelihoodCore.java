package dr.evomodel.newtreelikelihood;

import dr.app.beagle.evomodel.treelikelihood.AbstractTreeLikelihood;

import java.util.logging.Logger;

/**
 * @author Marc Suchard
 * @author Andrew Rambaut
 *
 */
public class ThreadedGPULikelihoodCore extends NativeLikelihoodCore {
    /**
     * Constructor
     *
     * @param instanceNumber identifies different TreeLikelihood instances
     * @param deviceNumber   specifies on which GPU device to run
     * @param stateCount     number of states
     * @param gpuInfo        complete information about GPUs
     */
    public ThreadedGPULikelihoodCore(int instanceNumber, int deviceNumber, int stateCount, GPUInfo gpuInfo) {

        StringBuffer sb = new StringBuffer();
        sb.append("Constructing GPU likelihood core: ");
        sb.append("Instance #" + (instanceNumber+1) + " on device #" + (deviceNumber+1) + "\n");
        sb.append("\tGPU Name: " + gpuInfo.getName(deviceNumber) + "\n");
        sb.append("If you publish results using this core, please cite Suchard and Rambaut (in preparation)\n");
        Logger.getLogger("dr.evomodel.treelikelihood").info(sb.toString());

        this.stateCount = stateCount;
        this.instanceNumber = instanceNumber;
        this.deviceNumber = deviceNumber;
    }

    public void finalize() throws Throwable {
        freeNativeMemory(instanceNumber);
    }

    public boolean canHandleDynamicRescaling() {
    	return true;
    }

    private int stateCount;
    private int instanceNumber;
    private int deviceNumber;

    public boolean canHandleTipPartials() {
        return true;
    }

    public boolean canHandleTipStates() {
        return (stateCount == 4);
    }

    public void initialize(int nodeCount, int stateTipCount, int patternCount, int matrixCount) {
        initialize(instanceNumber, deviceNumber, nodeCount, stateTipCount, patternCount, matrixCount, stateCount);
    }

    public void setTipPartials(int tipIndex, double[] partials) {
        setTipPartials(instanceNumber, tipIndex, partials);
    }

    public void setTipStates(int tipIndex, int[] states) {
        setTipStates(instanceNumber, tipIndex, states);
    }

    protected void updateRootFrequencies(double[] frequencies) {
        updateRootFrequencies(instanceNumber, frequencies);
    }

    protected void updateEigenDecomposition(double[][] eigenVectors, double[][] inverseEigenValues, double[] eigenValues) {
        updateEigenDecomposition(instanceNumber, eigenVectors, inverseEigenValues, eigenValues);
    }

    protected void updateCategoryRates(double[] rates) {
        updateCategoryRates(instanceNumber, rates);
    }

    protected void updateCategoryProportions(double[] proportions) {
        updateCategoryProportions(instanceNumber, proportions);
    }

    public void updateMatrices(int[] branchUpdateIndices, double[] branchLengths, int branchUpdateCount) {
        updateMatrices(instanceNumber, branchUpdateIndices, branchLengths, branchUpdateCount);
    }

    public void updatePartials(int[] operations, int[] dependencies, int operationCount, boolean rescale) {
        updatePartials(instanceNumber, operations, dependencies, operationCount, rescale);
    }

    public void calculateLogLikelihoods(int rootNodeIndex, double[] outLogLikelihoods) {
        calculateLogLikelihoods(instanceNumber, rootNodeIndex, outLogLikelihoods);
    }

    public void storeState() {
        storeState(instanceNumber);
    }

    public void restoreState() {
        restoreState(instanceNumber);
    }

    /**
     * GPU-specific loading *
     */

    protected static String GPU_LIBRARY_NAME = "ThreadedGPULikelihoodCore";

    protected static boolean isCompatible(GPUInfo gpuInfo, int[] configuration) {
        return true;
    }

    private static GPUInfo gpuInfo = null;

    private static int instanceCounter = 0;

    /**
     * GPU-specific native interface *
     */

    private static native GPUInfo getGPUInfo();

    private native void initialize(int instanceNumber, int deviceNumber, int nodeCount, int stateTipCount, int patternCount, int matrixCount, int stateCount);

    private native void freeNativeMemory(int instanceNumber);

    private native void setTipPartials(int instanceNumber, int tipIndex, double[] partials);

    private native void setTipStates(int instanceNumber, int tipIndex, int[] states);

    private native void updateRootFrequencies(int instanceNumber, double[] frequencies);

    private native void updateEigenDecomposition(int instanceNumber, double[][] eigenVectors, double[][] inverseEigenValues, double[] eigenValues);

    private native void updateCategoryRates(int instanceNumber, double[] rates);

    private native void updateCategoryProportions(int instanceNumber, double[] proportions);

    private native void updateMatrices(int instanceNumber, int[] branchUpdateIndices, double[] branchLengths, int branchUpdateCount);

    private native void updatePartials(int instanceNumber, int[] operations, int[] dependencies, int operationCount, boolean rescale);

    private native void calculateLogLikelihoods(int instanceNumber, int rootNodeIndex, double[] outLogLikelihoods);

    private native void storeState(int instanceNumber);

    private native void restoreState(int instanceNumber);

    /*
     * LikelihoodCoreLoader for LikelihoodCoreFactory *
     */

    public static class LikelihoodCoreLoader implements LikelihoodCoreFactory.LikelihoodCoreLoader {

        public String getLibraryName() {
            return GPU_LIBRARY_NAME;
        }

        public LikelihoodCore createLikelihoodCore(int[] configuration, AbstractTreeLikelihood treeLikelihood) {
            int stateCount = configuration[0];
            int paddedStateCount = stateCount;
            int deviceNumber = configuration[3];
            if (stateCount == 4)
                paddedStateCount = 4;
            else if (stateCount <= 16)
                paddedStateCount = 16;
            else if (stateCount <= 32)
                paddedStateCount = 32;
            else if (stateCount <= 64)
                paddedStateCount = 64;
            else if (stateCount >= 122 && stateCount <= 128)
            	paddedStateCount = 128;
        	else if (stateCount >= 183 && stateCount <= 192)
				paddedStateCount = 192;
            try {
                System.loadLibrary(getLibraryName() + "-" + paddedStateCount);
                if (gpuInfo == null) {
                    gpuInfo = ThreadedGPULikelihoodCore.getGPUInfo();
                    if (gpuInfo == null) // No GPU is present
                        return null;
                    Logger.getLogger("dr.evomodel.treelikelihood").info(gpuInfo.toString());
                    if (deviceNumber < 0 || deviceNumber >= gpuInfo.numberDevices) {
                        throw new RuntimeException("Cannot access GPU device #"+(deviceNumber+1));
                    }
                }
                if (!ThreadedGPULikelihoodCore.isCompatible(gpuInfo, configuration)) // GPU is not compatible
                    return null;
            } catch (UnsatisfiedLinkError e) { // No library present
                return null;
            }
            instanceCounter++;
            return new ThreadedGPULikelihoodCore(instanceCounter - 1, deviceNumber, stateCount, gpuInfo);
		}
	}
}
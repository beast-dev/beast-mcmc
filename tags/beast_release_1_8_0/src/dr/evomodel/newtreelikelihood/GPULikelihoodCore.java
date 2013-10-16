package dr.evomodel.newtreelikelihood;

import dr.app.beagle.evomodel.treelikelihood.AbstractTreeLikelihood;

import java.util.logging.Logger;

/**
 * @author Marc Suchard
 * @author Andrew Rambaut
 * @author Dat Hunyh
 *
 */

public class GPULikelihoodCore extends NativeLikelihoodCore {

	/**
	 * Constructor
	 *
	 * @param treeLikelihood reference back to ATL, should be able to remove at some point
	 * @param gpuInfo number of states
	 */
    public GPULikelihoodCore(int deviceNumber, int stateCount, AbstractTreeLikelihood treeLikelihood, GPUInfo gpuInfo) {

		StringBuffer sb = new StringBuffer();
		sb.append("Constructing GPU likelihood core:\n");
		sb.append("\tGPU Name: "+gpuInfo.getName(deviceNumber)+"\n");
//		sb.append("\tOther info here: "+gpuInfo.memorySize[deviceNumber]+"\n");
		sb.append("If you publish results using this core, please reference Suchard and Rambaut (in preparation)\n");
		Logger.getLogger("dr.evomodel.treelikelihood").info(sb.toString());
		this.treeLikelihood = treeLikelihood;
		this.stateCount = stateCount;
	}

    private int stateCount;

    private AbstractTreeLikelihood treeLikelihood;   // TODO Not needed if everything runs in one thread or I get GPU contexts working

	/** GPU-specific loading **/

	protected static String GPU_LIBRARY_NAME = "GPULikelihoodCore";

	protected static boolean isCompatible(GPUInfo gpuInfo, int[] configuration) {
		return true;
	}


	/** GPU-specific native interface **/

	private static native GPUInfo getGPUInfo();

    public boolean canHandleTipPartials() {
        return true;
    }

    public boolean canHandleTipStates() {
        return false;
    }

	/** Native interface overriding NativeLikelihoodCore **/

	public void initialize(int nodeCount, int stateTipCount, int patternCount, int matrixCount) {
//		System.err.println("stateCOunt = "+stateCount+" in java");
		initialize(nodeCount, stateTipCount, patternCount, matrixCount, stateCount);
	}

    public native void initialize(int nodeCount, int stateTipCount, int patternCount, int matrixCount, int stateCount);

    private native void freeNativeMemory();

    public native void setTipPartials(int tipIndex, double[] partials);

    public native void setTipStates(int tipIndex, int[] states);

	protected native void updateRootFrequencies(double[] frequencies);

	protected native void updateEigenDecomposition(double[][] eigenVectors, double[][] inverseEigenValues, double[] eigenValues);

	protected native void updateCategoryRates(double[] rates);

    protected native void updateCategoryProportions(double[] proportions);

    public native void updateMatrices(int[] branchUpdateIndices, double[] branchLengths, int branchUpdateCount);

    public native void updatePartials(int[] operations, int[] dependencies, int operationCount);

    public native void calculateLogLikelihoods(int rootNodeIndex, double[] outLogLikelihoods);

    public native void storeState();

    public native void restoreState();

    private static GPUInfo gpuInfo = null;

	public static class LikelihoodCoreLoader implements LikelihoodCoreFactory.LikelihoodCoreLoader {

		public String getLibraryName() { return GPU_LIBRARY_NAME; }

		public LikelihoodCore createLikelihoodCore(int[] configuration, AbstractTreeLikelihood treeLikelihood) {
			int stateCount = configuration[0];
			int paddedStateCount = stateCount;
			int deviceNumber = 0;
			if ( stateCount == 4 )
				paddedStateCount = 4;
			else if ( stateCount <= 16 )
				paddedStateCount = 16;
			else if (stateCount <= 32 )
				paddedStateCount = 32;
			else if (stateCount <= 64 )
				paddedStateCount = 64;
			try {
				System.loadLibrary(getLibraryName()+"-"+paddedStateCount);
				if (gpuInfo == null) {
					gpuInfo = GPULikelihoodCore.getGPUInfo();
					if (gpuInfo == null) // No GPU is present
						return null;
					Logger.getLogger("dr.evomodel.treelikelihood").info(gpuInfo.toString());
				}
				if (!GPULikelihoodCore.isCompatible(gpuInfo, configuration)) // GPU is not compatible
					return null;
			} catch (UnsatisfiedLinkError e) { // No library present
				return null;
			}
			return new GPULikelihoodCore(deviceNumber,stateCount, treeLikelihood, gpuInfo);
		}
	}
}
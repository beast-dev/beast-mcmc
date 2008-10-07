package dr.evomodel.treelikelihood;

import dr.evomodel.sitemodel.SiteModel;
import dr.math.matrixAlgebra.Vector;

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
		sb.append("\tGPU Name: \n");
		sb.append("\tOther info here\n");
		sb.append("If you publish results using this core, please reference Hunyh and Suchard (in preparation)\n");
		Logger.getLogger("dr.evomodel.treelikelihood").info(sb.toString());
	}

	public void finalize() throws Throwable {
		super.finalize();
		nativeFinalize();
	}

	private double[] tmpPartials = null;

	protected void calculateIntegratePartials(long ptrInPartials, double[] proportions, double[] outPartials) {

		if (tmpPartials == null)
			tmpPartials = new double[partialsSize];

		double[] inPartials = tmpPartials;

		getNativeMemoryArray(ptrInPartials, 0, tmpPartials, 0, partialsSize);

		int u = 0;
		int v = 0;
		for (int k = 0; k < patternCount; k++) {

			for (int i = 0; i < stateCount; i++) {

				outPartials[u] = inPartials[v] * proportions[0];
				u++;
				v++;
			}
		}


		for (int l = 1; l < matrixCount; l++) {
			u = 0;

			for (int k = 0; k < patternCount; k++) {

				for (int i = 0; i < stateCount; i++) {

					outPartials[u] += inPartials[v] * proportions[l];
					u++;
					v++;
				}
			}
		}

	}

	/* GPU memory handing functions */

	protected native long allocateNativeMemoryArray(int length);

	protected native long allocateNativeIntMemoryArray(int length);

	protected native void freeNativeMemoryArray(long nativePtr);

	protected native void setNativeMemoryArray(double[] fromArray, int fromOffset,
	                                           long toNativePtr, int toOffset, int length);

	protected native void setNativeMemoryArray(int[] fromArray, int fromOffset,
	                                           long toNativePtr, int toOffset, int length);

	protected native void getNativeMemoryArray(long fromNativePtr, int fromOffset,
	                                           double[] toArray, int toOffset, int length);

	protected native void getNativeMemoryArray(long fromNativePtr, int fromOffset,
	                                           int[] toArray, int toOffset, int length);

	protected native int getNativeRealSize();

	protected native int getNativeIntSize();

	protected native long getNativeRealDelta(long inPtr, int length);

	private native void nativeFinalize();

	/* GPU peeling functions */

	public native void nativeIntegratePartials(long ptrPartials, long ptrProportions,
	                                           long ptrOutPartials,
	                                           int patternCount, int matrixCount,
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
			System.loadLibrary("GPUMemoryLikelihoodCore");
			isNativeAvailable = true;
		} catch (UnsatisfiedLinkError e) {
		}
	}
}

package dr.evomodel.treelikelihood;

public class GPUGeneralLikelihoodCore extends AbstractLikelihoodCore {

	public GPUGeneralLikelihoodCore(int stateCount) {
		super(stateCount);
		String nameGPU = nativeCheckGPU();
		if (nameGPU == null) 
		    throw new RuntimeException("Missing GPU");
		System.err.println("Using GPU: "+nameGPU);
	}


	protected void calculateIntegratePartials(double[] inPartials,
	                                          double[] proportions, double[] outPartials) {
		nativeGPUIntegratePartials(inPartials, proportions, patternCount, matrixCount, outPartials, stateCount);
	}

	protected void calculatePartialsPartialsPruning(double[] partials1,
	                                                double[] matrices1, double[] partials2, double[] matrices2,
	                                                double[] partials3) {
		nativeGPUPartialsPartialsPruning(partials1, matrices1, partials2, matrices2, patternCount, matrixCount, partials3, stateCount);
	}

	protected void calculateStatesPartialsPruning(int[] states1,
	                                              double[] matrices1, double[] partials2, double[] matrices2,
	                                              double[] partials3) {
		nativeGPUStatesPartialsPruning(states1, matrices1, partials2, matrices2, patternCount, matrixCount, partials3, stateCount);
	}

	protected void calculateStatesStatesPruning(int[] states1,
	                                            double[] matrices1, int[] states2, double[] matrices2,
	                                            double[] partials3) {


		nativeGPUStatesStatesPruning(states1, matrices1, states2, matrices2, patternCount, matrixCount, partials3, stateCount);
	}

	protected void calculatePartialsPartialsPruning(double[] partials1,
	                                                double[] matrices1, double[] partials2, double[] matrices2,
	                                                double[] partials3, int[] matrixMap) {
		throw new RuntimeException("not implemented using matrixMap");
	}

	protected void calculateStatesStatesPruning(int[] states1,
	                                            double[] matrices1, int[] states2, double[] matrices2,
	                                            double[] partials3, int[] matrixMap) {
		throw new RuntimeException("not implemented using matrixMap");
	}

	protected void calculateStatesPartialsPruning(int[] states1,
	                                              double[] matrices1, double[] partials2, double[] matrices2,
	                                              double[] partials3, int[] matrixMap) {
		throw new RuntimeException("not implemented using matrixMap");
	}

	public void calculateLogLikelihoods(double[] partials,
	                                    double[] frequencies, double[] outLogLikelihoods) {
		int v = 0;
		for (int k = 0; k < patternCount; k++) {

			double sum = 0.0;
			for (int i = 0; i < stateCount; i++) {

				sum += frequencies[i] * partials[v];
				v++;
			}
			outLogLikelihoods[k] = Math.log(sum) + getLogScalingFactor(k);
		}
	}

	public native void nativeIntegratePartials(double[] partials, double[] proportions,
	                                           int patternCount, int matrixCount,
	                                           double[] outPartials,
	                                           int stateCount);

	public native void nativeGPUIntegratePartials(double[] partials, double[] proportions,
	                                           int patternCount, int matrixCount,
	                                           double[] outPartials,
	                                           int stateCount);

	protected native void nativePartialsPartialsPruning(double[] partials1, double[] matrices1,
	                                                    double[] partials2, double[] matrices2,
	                                                    int patternCount, int matrixCount,
	                                                    double[] partials3,
	                                                    int stateCount);

	protected native void nativeGPUPartialsPartialsPruning(double[] partials1, double[] matrices1,
	                                                    double[] partials2, double[] matrices2,
	                                                    int patternCount, int matrixCount,
	                                                    double[] partials3,
	                                                    int stateCount);

	protected native void nativeStatesPartialsPruning(int[] states1, double[] matrices1,
	                                                  double[] partials2, double[] matrices2,
	                                                  int patternCount, int matrixCount,
	                                                  double[] partials3,
	                                                  int stateCount);

	protected native void nativeGPUStatesPartialsPruning(int[] states1, double[] matrices1,
	                                                  double[] partials2, double[] matrices2,
	                                                  int patternCount, int matrixCount,
	                                                  double[] partials3,
	                                                  int stateCount);

	protected native void nativeStatesStatesPruning(int[] states1, double[] matrices1,
	                                                int[] states2, double[] matrices2,
	                                                int patternCount, int matrixCount,
	                                                double[] partials3,
	                                                int stateCount);

	protected native void nativeGPUStatesStatesPruning(int[] states1, double[] matrices1,
	                                                int[] states2, double[] matrices2,
	                                                int patternCount, int matrixCount,
	                                                double[] partials3,
	                                                int stateCount);


	protected native void setMode(int mode);

        protected native String nativeCheckGPU();

	public static boolean isAvailable() {
		return isNativeGPUAvailable;
	}

	private static boolean isNativeGPUAvailable = false;

	static {
		try {
			System.loadLibrary("GPUGeneralLikelihoodCore");
			isNativeGPUAvailable = true;
		} catch (UnsatisfiedLinkError e) {
//            System.err.println("Using Java general likelihood core " + e.toString());
		}
	}
}
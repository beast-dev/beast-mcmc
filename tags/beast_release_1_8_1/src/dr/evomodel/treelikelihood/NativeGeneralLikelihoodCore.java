package dr.evomodel.treelikelihood;

public class NativeGeneralLikelihoodCore extends AbstractLikelihoodCore {

	public NativeGeneralLikelihoodCore(int stateCount) {
		super(stateCount);
	}


	protected void calculateIntegratePartials(double[] inPartials,
	                                          double[] proportions, double[] outPartials) {
		nativeIntegratePartials(inPartials, proportions, patternCount, matrixCount, outPartials, stateCount);
	}

	protected void calculatePartialsPartialsPruning(double[] partials1,
	                                                double[] matrices1, double[] partials2, double[] matrices2,
	                                                double[] partials3) {
		nativePartialsPartialsPruning(partials1, matrices1, partials2, matrices2, patternCount, matrixCount, partials3, stateCount);
	}

	protected void calculateStatesPartialsPruning(int[] states1,
	                                              double[] matrices1, double[] partials2, double[] matrices2,
	                                              double[] partials3) {
		nativeStatesPartialsPruning(states1, matrices1, partials2, matrices2, patternCount, matrixCount, partials3, stateCount);
	}

	protected void calculateStatesStatesPruning(int[] states1,
	                                            double[] matrices1, int[] states2, double[] matrices2,
	                                            double[] partials3) {


		nativeStatesStatesPruning(states1, matrices1, states2, matrices2, patternCount, matrixCount, partials3, stateCount);
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

	protected native void nativePartialsPartialsPruning(double[] partials1, double[] matrices1,
	                                                    double[] partials2, double[] matrices2,
	                                                    int patternCount, int matrixCount,
	                                                    double[] partials3,
	                                                    int stateCount);

	protected native void nativeStatesPartialsPruning(int[] states1, double[] matrices1,
	                                                  double[] partials2, double[] matrices2,
	                                                  int patternCount, int matrixCount,
	                                                  double[] partials3,
	                                                  int stateCount);

	protected native void nativeStatesStatesPruning(int[] states1, double[] matrices1,
	                                                int[] states2, double[] matrices2,
	                                                int patternCount, int matrixCount,
	                                                double[] partials3,
	                                                int stateCount);

	protected native void setMode(int mode);

	public static boolean isAvailable() {
		return isNativeAvailable;
	}

	private static boolean isNativeAvailable = false;

	static {
		try {
			System.loadLibrary("GeneralLikelihoodCore");
			isNativeAvailable = true;
		} catch (UnsatisfiedLinkError e) {
			System.err.println("Unable to load general core");
		}
	}
}
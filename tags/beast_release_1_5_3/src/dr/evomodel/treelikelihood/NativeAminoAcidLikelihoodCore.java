package dr.evomodel.treelikelihood;

public class NativeAminoAcidLikelihoodCore extends AbstractLikelihoodCore{

	public NativeAminoAcidLikelihoodCore() {
		super(20);
	}


	protected void calculateIntegratePartials(double[] inPartials,
			double[] proportions, double[] outPartials) {
		nativeIntegratePartials(inPartials, proportions, patternCount, matrixCount, outPartials);
	}

	protected void calculatePartialsPartialsPruning(double[] partials1,
			double[] matrices1, double[] partials2, double[] matrices2,
			double[] partials3) {
		nativePartialsPartialsPruning(partials1, matrices1, partials2, matrices2, patternCount, matrixCount, partials3);
	}

	protected void calculateStatesPartialsPruning(int[] states1,
			double[] matrices1, double[] partials2, double[] matrices2,
			double[] partials3) {
		nativeStatesPartialsPruning(states1, matrices1, partials2, matrices2, patternCount, matrixCount, partials3);
	}

	protected void calculateStatesStatesPruning(int[] states1,
			double[] matrices1, int[] states2, double[] matrices2,
			double[] partials3) {


		nativeStatesStatesPruning(states1, matrices1, states2, matrices2, patternCount, matrixCount, partials3);
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

			double sum = frequencies[0] * partials[v];	v++;
			sum += frequencies[1] * partials[v];	v++;
			sum += frequencies[2] * partials[v];	v++;
			sum += frequencies[3] * partials[v];	v++;

			sum += frequencies[4] * partials[v];	v++;
			sum += frequencies[5] * partials[v];	v++;
			sum += frequencies[6] * partials[v];	v++;
			sum += frequencies[7] * partials[v];	v++;

			sum += frequencies[8] * partials[v];	v++;
			sum += frequencies[9] * partials[v];	v++;
			sum += frequencies[10] * partials[v];	v++;
			sum += frequencies[11] * partials[v];	v++;

			sum += frequencies[12] * partials[v];	v++;
			sum += frequencies[13] * partials[v];	v++;
			sum += frequencies[14] * partials[v];	v++;
			sum += frequencies[15] * partials[v];	v++;

			sum += frequencies[16] * partials[v];	v++;
			sum += frequencies[17] * partials[v];	v++;
			sum += frequencies[18] * partials[v];	v++;
			sum += frequencies[19] * partials[v];	v++;
            outLogLikelihoods[k] = Math.log(sum) + getLogScalingFactor(k);
		}
	}

	public native void nativeIntegratePartials(double[] partials, double[] proportions,
											   int patternCount, int matrixCount,
											   double[] outPartials);

	protected native void nativePartialsPartialsPruning(double[] partials1, double[] matrices1,
														double[] partials2, double[] matrices2,
														int patternCount, int matrixCount,
														double[] partials3);

	protected native void nativeStatesPartialsPruning(int[] states1, double[] matrices1,
													  double[] partials2, double[] matrices2,
													  int patternCount, int matrixCount,
													  double[] partials3);

	protected native void nativeStatesStatesPruning(int[] states1, double[] matrices1,
													int[] states2, double[] matrices2,
													int patternCount, int matrixCount,
													double[] partials3);

	public static boolean isAvailable(){
		return isNativeAvailable;
	}

	private static boolean isNativeAvailable = false;

	static {
        try {
            System.loadLibrary("AminoAcidLikelihoodCore");
            isNativeAvailable = true;
        } catch (UnsatisfiedLinkError e) {
            System.err.println("Using Java AminoAcid likelihood core " + e.toString());
            System.err.println("Looking for AminoAcidLikelihoodCore in " + System.getProperty("java.library.path"));
        }
    }
}

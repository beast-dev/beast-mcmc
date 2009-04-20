package dr.evomodel.branchratemodel;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.inference.model.Parameter;

/**
 * @author Marc Suchard
 */
public class ContinuousEpochBranchRateModel extends RateEpochBranchRateModel {

	/**
	 * The constructor. For an N-epoch model, there should be N rate paramters and N-1 transition times.
	 *
	 * @param timeParameters an array of transition time parameters
	 * @param rateParameters an array of rate parameters
	 */
	public ContinuousEpochBranchRateModel(Parameter[] timeParameters, Parameter[] rateParameters, Parameter rootHeight) {
		super(timeParameters, rateParameters);
		this.rootHeight = rootHeight;
		addParameter(rootHeight);
		normalizationKnown = false;
	}

	protected void handleParameterChangedEvent(Parameter parameter, int index) {
		fireModelChanged();
		normalizationKnown = false;
	}

	private void normalize() {
		normalization = 0.0;
		double nextTime;
		double thisRate;
		double startTime = 0.0;
		for (int j = 0; j < timeParameters.length; j++) {
			nextTime = timeParameters[j].getParameterValue(0);
			thisRate = rateParameters[j].getParameterValue(0);
			normalization += (nextTime - startTime) * thisRate;
			startTime = nextTime;
		}
		nextTime = rootHeight.getParameterValue(0);
		thisRate = rateParameters[timeParameters.length].getParameterValue(0);
		normalization += (nextTime - startTime) * thisRate;
	}

	protected void storeState() {
		savedNormalization = normalization;
	}

	protected void restoreState() {
		normalization = savedNormalization;
	}

	public double getBranchRate(Tree tree, NodeRef node) {

		if (!normalizationKnown)
			normalize();

		NodeRef parent = tree.getParent(node);

		if (parent != null) {
			double height0 = tree.getNodeHeight(node);
			double height1 = tree.getNodeHeight(parent);
			int i = 0;

			double rate = 0.0;
			double lastHeight = height0;

			// First find the epoch which contains the node height
			while (i < timeParameters.length && height0 > timeParameters[i].getParameterValue(0)) {
				i++;
			}

			// Now walk up the branch until we reach the last epoch or the height of the parent
			while (i < timeParameters.length && height1 > timeParameters[i].getParameterValue(0)) {
				// add the rate for that epoch multiplied by the time spent at that rate
				rate += rateParameters[i].getParameterValue(0) * (timeParameters[i].getParameterValue(0) - lastHeight);
				lastHeight = timeParameters[i].getParameterValue(0);
				i++;
			}

			// Add that last rate segment
			rate += rateParameters[i].getParameterValue(0) * (height1 - lastHeight);

			return rate / normalization / (height1 - height0);
		}
		throw new IllegalArgumentException("root node doesn't have a rate!");
	}

	private Parameter rootHeight;
	private double normalization;
	private double savedNormalization;
	private boolean normalizationKnown = false;

}

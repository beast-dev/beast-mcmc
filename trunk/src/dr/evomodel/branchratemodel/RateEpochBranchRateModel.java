package dr.evomodel.branchratemodel;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.inference.model.*;
import dr.xml.*;

import java.util.*;
import java.util.logging.Logger;

/**
 * Implements a model where time is broken into 'epochs' each with a different but
 * constant rate.
 * @author Andrew Rambaut
 *
 * @version $Id$
 */
public class RateEpochBranchRateModel extends AbstractModel implements BranchRateModel  {

	public static final String RATE_EPOCH_BRANCH_RATES = "rateEpochBranchRates";
	public static final String RATE = "rate";
	public static final String EPOCH = "epoch";
	public static final String TRANSITION_TIME = "transitionTime";

	private final double[] epochTransitionTimes;
	private final Parameter[] rateParameters;

	/**
	 * The constructor. For an N-epoch model, there should be N rate paramters and N-1 transition times.
	 * @param epochTransitionTimes an array of transition times
	 * @param rateParameters an array of rate parameters
	 */
	public RateEpochBranchRateModel(double[] epochTransitionTimes, Parameter[] rateParameters) {

		super(RATE_EPOCH_BRANCH_RATES);

		this.epochTransitionTimes = epochTransitionTimes;
		this.rateParameters = rateParameters;
		for (Parameter parameter : rateParameters) {
			addParameter(parameter);
		}
	}

	public void handleModelChangedEvent(Model model, Object object, int index) {
		// nothing to do
	}

	protected void handleParameterChangedEvent(Parameter parameter, int index) {
		fireModelChanged();
	}

	protected void storeState() {
		// nothing to do
	}

	protected void restoreState() {
		// nothing to do
	}

	protected void acceptState() {
		// nothing to do
	}

	public double getBranchRate(Tree tree, NodeRef node) {

		NodeRef parent = tree.getParent(node);

		if (parent != null) {
			double height0 = tree.getNodeHeight(node);
			double height1 = tree.getNodeHeight(parent);
			int i = 0;

			double rate = 0.0;
			double lastHeight = height0;

			// First find the epoch which contains the node height
			while (i < epochTransitionTimes.length && height0 > epochTransitionTimes[i]) {
				i++;
			}

			// Now walk up the branch until we reach the last epoch or the height of the parent
			while (i < epochTransitionTimes.length && height1 > epochTransitionTimes[i]) {
				// add the rate for that epoch multiplied by the time spent at that rate
				rate += rateParameters[i].getParameterValue(0) * (epochTransitionTimes[i] - lastHeight);
				lastHeight = epochTransitionTimes[i];
				i++;
			}

			// Add that last rate segment
			rate += rateParameters[i].getParameterValue(0) * (height1 - lastHeight);

			// normalize the rate for the branch length
			rate = rate / (height1 - height0);

			return rate;
		}
		throw new IllegalArgumentException("root node doesn't have a rate!");
	}

	public String getBranchAttributeLabel() {
		return "rate";
	}

	public String getAttributeForBranch(Tree tree, NodeRef node) {
		return Double.toString(getBranchRate(tree, node));
	}

	public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

		public String getParserName() { return RATE_EPOCH_BRANCH_RATES; }

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {

			Logger.getLogger("dr.evomodel").info("Using multi-epoch rate model.");

			List<Epoch> epochs = new ArrayList<Epoch>();

			for (int i = 0; i < xo.getChildCount(); i++) {
				XMLObject xoc = (XMLObject)xo.getChild(i);
				if (xoc.getName().equals(EPOCH)) {
					double t = xoc.getDoubleAttribute(TRANSITION_TIME);
					Parameter p = (Parameter)xoc.getChild(Parameter.class);
					epochs.add(new Epoch(t, p));
				}
			}

			Parameter ancestralRateParameter = (Parameter)xo.getSocketChild(RATE);

			Collections.sort(epochs);
			double[] epochTransitionTimes = new double[epochs.size()];
			Parameter[] rateParameters = new Parameter[epochs.size() + 1];

			int i = 0;
			for (Epoch epoch : epochs) {
				epochTransitionTimes[i] = epoch.transitionTime;
				rateParameters[i] = epoch.parameter;
				i++;
			}
			rateParameters[i] = ancestralRateParameter;

			return new RateEpochBranchRateModel(epochTransitionTimes, rateParameters);
		}

		class Epoch implements Comparable {

			private final double transitionTime;
			private final Parameter parameter;

			public Epoch(double transitionTime, Parameter parameter) {
				this.transitionTime = transitionTime;
				this.parameter = parameter;
			}

			public int compareTo(Object o) {
				return Double.compare(transitionTime, ((Epoch)o).transitionTime);
			}

		}
		//************************************************************************
		// AbstractXMLObjectParser implementation
		//************************************************************************

		public String getParserDescription() {
			return
					"This element provides a multiple epoch molecular clock model. " +
							"All branches (or portions of them) have the same rate of molecular evolution within a given epoch.";
		}

		public Class getReturnType() { return RateEpochBranchRateModel.class; }

		public XMLSyntaxRule[] getSyntaxRules() { return rules; }

		private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
				new ElementRule(EPOCH,
						new XMLSyntaxRule[] {
								AttributeRule.newDoubleRule(TRANSITION_TIME, false, "The time of transition between this epoch and the previous one"),
								new ElementRule(Parameter.class, "The evolutionary rate parameter for this epoch"),
						}, "An epoch that lasts until transitionTime",
						1, Integer.MAX_VALUE
				),
				new ElementRule(RATE, Parameter.class, "The ancestral molecular evolutionary rate parameter", false)
		};
	};

}

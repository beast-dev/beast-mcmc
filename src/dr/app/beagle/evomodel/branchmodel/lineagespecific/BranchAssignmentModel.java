package dr.app.beagle.evomodel.branchmodel.lineagespecific;

import dr.app.beagle.evomodel.branchmodel.BranchModel;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.AbstractBranchRateModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Model;
import dr.inference.model.Variable;
import dr.inference.model.Variable.ChangeType;
import dr.math.GammaFunction;
import dr.math.distributions.Distribution;

public class BranchAssignmentModel extends AbstractBranchRateModel {

	private Distribution discreteDistribution;
	private int categoryCount;

	private TreeModel treeModel;
	private BranchModel branchSpecificModel;

	public BranchAssignmentModel(BranchModel branchSpecificModel,
			TreeModel treeModel, Distribution discreteDistribution, int categoryCount) {

		super("");
		
		this.branchSpecificModel = branchSpecificModel;
		this.treeModel = treeModel;
		this.discreteDistribution = discreteDistribution;
		this.categoryCount = categoryCount;

	}// END: Constructor

	private int[] getBranchAssignmentCounts() {

		int[] branchAssignments = new int[categoryCount];

		for (NodeRef branch : treeModel.getNodes()) {

			int branchParameterIndex = ((BranchSpecific) branchSpecificModel)
					.getBranchModelMapping(branch).getOrder()[0];
			branchAssignments[branchParameterIndex]++;
		}

		return branchAssignments;
	}// END: getBranchAssignmentCounts

	/*
	 * Multinomial Distribution Likelihood
	 */
	public double getLogLikelihood() {

		double logLike = 0.0;
		int[] branchAssignments = getBranchAssignmentCounts();
		double[] probs = ((StickBreakingProcessPrior) discreteDistribution)
				.getProbabilities();

		logLike += GammaFunction.lnGamma(categoryCount + 1);
		for (int i = 0; i < categoryCount; i++) {

			logLike += ((probs[i] * Math.log(branchAssignments[i] + 1)) - (GammaFunction
					.lnGamma(branchAssignments[i] + 1)));

		}// END: i loop

		return logLike;
	}// END: getLogLikelihood

	@Override
	public double getBranchRate(Tree tree, NodeRef node) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	protected void handleModelChangedEvent(Model model, Object object, int index) {
	}

	@Override
	protected void handleVariableChangedEvent(Variable variable, int index,
			ChangeType type) {
		
	}

	@Override
	protected void storeState() {
		
	}

	@Override
	protected void restoreState() {
		
	}

	@Override
	protected void acceptState() {
	}

}// END: class

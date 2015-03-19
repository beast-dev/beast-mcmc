package dr.app.beagle.evomodel.branchmodel.lineagespecific;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;

import dr.evolution.tree.NodeRef;
import dr.evomodel.branchratemodel.CountableBranchCategoryProvider;
import dr.evomodel.newtreelikelihood.TreeLikelihood;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.CompoundLikelihood;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;

public class BeagleBranchLikelihoods {

	TreeModel treeModel;
	List<Likelihood> treeLikelihoods;
	Parameter zParameter;

	// for discrete categories
	private CountableBranchCategoryProvider categoriesProvider;

	public BeagleBranchLikelihoods(TreeModel treeModel,
			List<Likelihood> treeLikelihoods, Parameter zParameter) {

		// super(treeLikelihoods);

		this.treeModel = treeModel;
		this.treeLikelihoods = treeLikelihoods;
		this.zParameter = zParameter;

		this.categoriesProvider = new CountableBranchCategoryProvider.IndependentBranchCategoryModel(
				treeModel, zParameter);

	}// END: Constructor

	public List<Likelihood> getBranchLikelihoods() {

		List<Likelihood> branchLikelihoods = new ArrayList<Likelihood>();

		for (NodeRef branch : treeModel.getNodes()) {

			if (!treeModel.isRoot(branch)) {

				int branchCategory = categoriesProvider.getBranchCategory(
						treeModel, branch);
				int zIndex = (int) zParameter.getParameterValue(branchCategory);

				branchLikelihoods.add(treeLikelihoods.get(zIndex));

			}
		}// END: branch loop

		return branchLikelihoods;
	}// END: getBranchLikelihoods

}// END: class

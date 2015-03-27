package dr.app.beagle.evomodel.branchmodel.lineagespecific;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;

import dr.evolution.tree.NodeRef;
import dr.evomodel.branchratemodel.CountableBranchCategoryProvider;
import dr.evomodel.newtreelikelihood.TreeLikelihood;
import dr.evomodel.tree.TreeModel;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.NumberColumn;
import dr.inference.model.CompoundLikelihood;
import dr.inference.model.Likelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;

public class BeagleBranchLikelihood implements Likelihood {

	TreeModel treeModel;
	List<Likelihood> uniqueLikelihoods;
	
	// size of z:
	List<Likelihood> branchLikelihoods;
	Parameter zParameter;

	private String id = null;

	// for discrete categories
	private CountableBranchCategoryProvider categoriesProvider;

	public BeagleBranchLikelihood(TreeModel treeModel,
			List<Likelihood> likelihoods, Parameter zParameter) {

		// super(treeLikelihoods);

		this.treeModel = treeModel;
		this.uniqueLikelihoods = likelihoods;
		this.zParameter = zParameter;

		this.categoriesProvider = new CountableBranchCategoryProvider.IndependentBranchCategoryModel(
				treeModel, zParameter);

		this.branchLikelihoods = getBranchLikelihoods();
		
	}// END: Constructor

	public List<Likelihood> getBranchLikelihoods() {

		List<Likelihood>  loglikes = new ArrayList<Likelihood>();

		for (NodeRef branch : treeModel.getNodes()) {

			if (!treeModel.isRoot(branch)) {

				int branchCategory = categoriesProvider.getBranchCategory(
						treeModel, branch);
				int zIndex = (int) zParameter.getParameterValue(branchCategory);

				Likelihood branchLikelihood = uniqueLikelihoods.get(zIndex);

				// branchLikelihoods.add(new Holder(branchLikelihood).value);
				loglikes.add(branchLikelihood);

			}
		}// END: branch loop

		return loglikes;
	}// END: getBranchLikelihoods

	// ///////////////////////
	// ---PRIVATE METHODS---//
	// ///////////////////////

	//

	// //////////////////////
	// ---PUBLIC METHODS---//
	// //////////////////////

	@Override
	public double getLogLikelihood() {
		
		double loglike = 0;
		for(Likelihood like : getBranchLikelihoods()) {
			
			loglike += like.getLogLikelihood();
			
		}//END: loglikes loop
		
		return loglike;
	}//END: getLogLikelihood

    public int getLikelihoodCount() {
        return branchLikelihoods.size();
    }

    public final Likelihood getLikelihood(int i) {
        return branchLikelihoods.get(i);
    }

    public List<Likelihood> getLikelihoods() {
        return branchLikelihoods;
    }
	
	@Override
	public LogColumn[] getColumns() {
		return new dr.inference.loggers.LogColumn[] { new LikelihoodColumn(
				getId() == null ? "likelihood" : getId()) };
	}

	@Override
	public String getId() {
		return this.id;
	}

	@Override
	public void setId(String id) {
		this.id = id;
	}

	@Override
	public Model getModel() {
		// TODO Auto-generated method stub
		return null;
	}


	@Override
	public void makeDirty() {
		// TODO Auto-generated method stub

	}

	@Override
	public String prettyName() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public boolean isUsed() {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public void setUsed() {
		// TODO Auto-generated method stub

	}

	@Override
	public boolean evaluateEarly() {
		// TODO Auto-generated method stub
		return false;
	}

	// ///////////////////////
	// ---PRIVATE CLASSES---//
	// ///////////////////////

	private class LikelihoodColumn extends NumberColumn {

		public LikelihoodColumn(String label) {
			super(label);
		}// END: Constructor

		public double getDoubleValue() {
			return getLogLikelihood();
		}

	}// END: LikelihoodColumn class

	// class Holder {
	//
	// public Likelihood value;
	//
	// public Holder(Likelihood initial) {
	// this.value = initial;
	// }// END: Holder
	//
	// }// END: Holder class

}// END: class

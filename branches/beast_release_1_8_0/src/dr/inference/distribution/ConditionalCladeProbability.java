package dr.inference.distribution;

import dr.evolution.tree.SimpleTree;
import dr.evomodel.tree.ConditionalCladeFrequency;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Likelihood;

public class ConditionalCladeProbability extends Likelihood.Abstract {
	
	private ConditionalCladeFrequency ccf;
	private TreeModel treeModel;
	
	public ConditionalCladeProbability(ConditionalCladeFrequency ccf, TreeModel treeModel) {
		super(null);
		this.ccf = ccf;
		this.treeModel = treeModel;
	}

	@Override
	protected double calculateLogLikelihood() {
		
		SimpleTree simTree = new SimpleTree(treeModel);
		
		//System.err.println("tree: " + simTree);
		//System.err.println(ccf.getTreeProbability(simTree));
		
		return ccf.getTreeProbability(simTree);
		
	}
	
	/**
     * Overridden to always return false.
     */
    protected boolean getLikelihoodKnown() {
        return false;
    }

}

package dr.evomodel.graphlikelihood;

import dr.evolution.alignment.PatternList;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.graph.GraphModel;
import dr.evomodel.graph.PartitionModel;
import dr.evomodel.graph.PartitionModel.PartitionChangedEvent;
import dr.evomodel.sitemodel.SiteModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treelikelihood.TipPartialsModel;
import dr.evomodel.treelikelihood.TreeLikelihood;
import dr.inference.model.Model;

public class GraphLikelihood extends TreeLikelihood {

    public GraphLikelihood(PatternList patternList,
            GraphModel graphModel,
            SiteModel siteModel,
            BranchRateModel branchRateModel,
            TipPartialsModel tipPartialsModel,
            boolean useAmbiguities,
            boolean allowMissingTaxa,
            boolean storePartials,
            boolean forceJavaCore,
            boolean forceRescaling) 
    {
    	super(patternList,graphModel,siteModel,branchRateModel,
    			tipPartialsModel,useAmbiguities,allowMissingTaxa,
    			storePartials,forceJavaCore,forceRescaling);
    	updatePattern = new boolean[patternList.getPatternCount()];
    }

    protected void growNodeStorage(){
    	boolean[] tmp = new boolean[updateNode.length*2];
    	System.arraycopy(updateNode, 0, tmp, 0, updateNode.length);
    	for(int i=updateNode.length; i<tmp.length; i++)
    		tmp[i]=false;
    	updateNode = tmp;
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
    	if(model instanceof GraphModel && 
    			object instanceof TreeModel.TreeChangedEvent &&
    			((TreeModel.TreeChangedEvent)object).getNode().getNumber() >= nodeCount)
    	{
    		growNodeStorage(); // need more node storage!
    	}
    	if(model instanceof PartitionModel){
    		// the partitions changed, do something here
    		// just mark all affected site patterns as dirty
    		PartitionChangedEvent pce = (PartitionChangedEvent)object;
    		if(pce.hasNewSection())
    		{
    			// mark only the new section as dirty
    			for(int i=pce.getNewSectionLeft(); i<pce.getNewSectionRight(); i++)
    			{
    				updatePattern[i]=true;
    			}
    		}
    	}else{
    		// otherwise the TreeLikelihood can handle it
            super.handleModelChangedEvent(model, object, index);
    	}
    }

    /**
     * partition-aware traversal
     */
    protected boolean traverse(Tree tree, NodeRef node) 
    {
    }
}

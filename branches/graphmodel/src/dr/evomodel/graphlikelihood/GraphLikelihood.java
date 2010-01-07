package dr.evomodel.graphlikelihood;

import java.util.logging.Logger;

import jebl.gui.trees.treeviewer.treelayouts.AbstractTreeLayout;
import dr.evolution.alignment.SimpleSiteList;
import dr.evolution.alignment.SiteList;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.TaxonList;
import dr.evomodel.graph.GraphModel;
import dr.evomodel.graph.PartitionModel;
import dr.evomodel.graph.PartitionModel.PartitionChangedEvent;
import dr.evomodel.sitemodel.SiteModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treelikelihood.AbstractTreeLikelihood;
import dr.evomodel.treelikelihood.GeneralLikelihoodCore;
import dr.evomodel.treelikelihood.LikelihoodCore;
import dr.evomodel.treelikelihood.TipPartialsModel;
import dr.evomodel.treelikelihood.TreeLikelihood;
import dr.inference.model.Model;

public class GraphLikelihood extends AbstractTreeLikelihood {

	public static final String GRAPH_LIKELIHOOD = "graphLikelihood";
    public static final String USE_AMBIGUITIES = "useAmbiguities";
    public static final String ALLOW_MISSING_TAXA = "allowMissingTaxa";
    public static final String STORE_PARTIALS = "storePartials";
    public static final String SCALING_FACTOR = "scalingFactor";
    public static final String SCALING_THRESHOLD = "scalingThreshold";
    public static final String FORCE_JAVA_CORE = "forceJavaCore";
    public static final String FORCE_RESCALING = "forceRescaling";
	SiteList concatSiteList;
	
    public GraphLikelihood(
            GraphModel graphModel,
            PartitionModel partitionModel,
            TipPartialsModel tipPartialsModel,
            boolean useAmbiguities,
            boolean allowMissingTaxa,
            boolean forceRescaling) 
    {
        super(GRAPH_LIKELIHOOD, createConcatenatedSiteList(partitionModel), graphModel);

        try{

        this.concatSiteList = (SiteList)this.patternList;

        // BEGIN TreeLikelihood LIKELIHOODCORE INIT CODE
        this.likelihoodCore = new GeneralLikelihoodCore(concatSiteList.getSiteCount());
        String coreName = "Java general";
        // END TreeLikelihood LIKELIHOODCORE INIT CODE

        SiteModel siteModel = (SiteModel)partitionModel.getModelsOnPartition(partitionModel.getPartition(0)).get(0);
        this.integrateAcrossCategories = siteModel.integrateAcrossCategories();
        this.categoryCount = siteModel.getCategoryCount();

    	updatePattern = new boolean[patternList.getPatternCount()];

    	
        // BEGIN TreeLikelihood LOGGER CODE
        final Logger logger = Logger.getLogger("dr.evomodel");
        {
            final String id = getId();
            logger.info("TreeLikelihood(" + ((id != null) ? id : treeModel.getId()) + ") using " + coreName + " likelihood core");

            logger.info("  " + (useAmbiguities ? "Using" : "Ignoring") + " ambiguities in tree likelihood.");
            logger.info("  With " + patternList.getPatternCount() + " unique site patterns.");
        }
        // END TreeLikelihood LOGGER CODE

        likelihoodCore.initialize(nodeCount, patternCount, categoryCount, integrateAcrossCategories);
        probabilities = new double[stateCount * stateCount];

    
        // BEGIN TreeLikelihood TIPPARTIALS CODE (VERBATIM)
        int extNodeCount = treeModel.getExternalNodeCount();
        int intNodeCount = treeModel.getInternalNodeCount();

        if (tipPartialsModel != null) {
            tipPartialsModel.setTree(treeModel);

            tipPartials = new double[patternCount * stateCount];

            for (int i = 0; i < extNodeCount; i++) {
                // Find the id of tip i in the patternList
                String id = treeModel.getTaxonId(i);
                int index = patternList.getTaxonIndex(id);

                if (index == -1) {
                    throw new TaxonList.MissingTaxonException("Taxon, " + id + ", in tree, " + treeModel.getId() +
                            ", is not found in patternList, " + patternList.getId());
                }

                tipPartialsModel.setStates(patternList, index, i, id);
                likelihoodCore.createNodePartials(i);
            }

            addModel(tipPartialsModel);
            //useAmbiguities = true;
        } else {
            for (int i = 0; i < extNodeCount; i++) {
                // Find the id of tip i in the patternList
                String id = treeModel.getTaxonId(i);
                int index = patternList.getTaxonIndex(id);

                if (index == -1) {
                    if (!allowMissingTaxa) {
                        throw new TaxonList.MissingTaxonException("Taxon, " + id + ", in tree, " + treeModel.getId() +
                                ", is not found in patternList, " + patternList.getId());
                    }
                    if (useAmbiguities) {
                        setMissingPartials(likelihoodCore, i);
                    } else {
                        setMissingStates(likelihoodCore, i);
                    }
                } else {
                    if (useAmbiguities) {
                        setPartials(likelihoodCore, patternList, categoryCount, index, i);
                    } else {
                        setStates(likelihoodCore, patternList, index, i);
                    }
                }
            }
        }
        for (int i = 0; i < intNodeCount; i++) {
            likelihoodCore.createNodePartials(extNodeCount + i);
        }

        if (forceRescaling) {
            likelihoodCore.setUseScaling(true);
            logger.info("  Forcing use of partials rescaling.");
        }
        // END TreeLikelihood TIPPARTIALS CODE (VERBATIM)

    
	    } catch (TaxonList.MissingTaxonException mte) {
	        throw new RuntimeException(mte.toString());
	    }

	}
    
    protected static SiteList createConcatenatedSiteList(PartitionModel pm){
    	SimpleSiteList ssl = new SimpleSiteList(pm.getSiteList(0).getDataType());
    	for(int i=0; i<pm.getSiteListCount(); i++){
    		SiteList sl = pm.getSiteList(i);
    		for(int j=0; j<sl.getSiteCount(); j++){
        		ssl.addPattern(sl.getSitePattern(j));
    		}
    	}
    	return ssl;
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

    protected void storeState() {
        likelihoodCore.storeState();
        super.storeState();
    }
    protected void restoreState() {
        likelihoodCore.restoreState();
        super.restoreState();
    }

    /**
     * Calculate the log likelihood of the current state.
     * Ripped nearly verbatim from treeLikelihood
     * @return the log likelihood.
     */
    protected double calculateLogLikelihood() {

        if (patternLogLikelihoods == null) {
            patternLogLikelihoods = new double[patternCount];
        }



        final NodeRef root = treeModel.getRoot();
        traverse(treeModel, root);

        double logL = 0.0;
        for (int i = 0; i < patternCount; i++) {
            logL += (patternLogLikelihoods[i]) * patternWeights[i];
        }

        if (logL == Double.NEGATIVE_INFINITY) {
            Logger.getLogger("dr.evomodel").info("TreeLikelihood, " + this.getId() + ", turning on partial likelihood scaling to avoid precision loss");

            // We probably had an underflow... turn on scaling
            likelihoodCore.setUseScaling(true);

            // and try again...
            updateAllNodes();
            updateAllPatterns();
            traverse(treeModel, root);

            logL = 0.0;
            for (int i = 0; i < patternCount; i++) {
                logL += (patternLogLikelihoods[i]) * patternWeights[i];
            }
        }

        //********************************************************************
        // after traverse all nodes and patterns have been updated --
        //so change flags to reflect this.
        for (int i = 0; i < nodeCount; i++) {
            updateNode[i] = false;
        }
        //********************************************************************

        return logL;
    }
    
    
    /**
     * partition-aware traversal
     */
    protected boolean traverse(Tree tree, NodeRef node) 
    {
    	return true;
    }


    protected GeneralLikelihoodCore likelihoodCore;
    protected int categoryCount;
    protected final boolean integrateAcrossCategories;
    protected double[] probabilities;
    protected double[] tipPartials;
    protected double[] patternLogLikelihoods = null;
    protected int[] siteCategories = null;
}

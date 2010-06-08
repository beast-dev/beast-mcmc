package dr.evomodel.graphlikelihood;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.logging.Logger;

import jebl.gui.trees.treeviewer.treelayouts.AbstractTreeLayout;
import dr.evolution.alignment.SimpleSiteList;
import dr.evolution.alignment.SiteList;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.TaxonList;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.graph.GraphModel;
import dr.evomodel.graph.Partition;
import dr.evomodel.graph.PartitionModel;
import dr.evomodel.graph.PartitionModel.PartitionChangedEvent;
import dr.evomodel.sitemodel.SiteModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.tree.TreeModel.TreeChangedEvent;
import dr.evomodel.treelikelihood.AbstractTreeLikelihood;
import dr.evomodel.treelikelihood.GeneralLikelihoodCore;
import dr.evomodel.treelikelihood.LikelihoodCore;
import dr.evomodel.treelikelihood.TipPartialsModel;
import dr.evomodel.treelikelihood.TreeLikelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;

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

	/*
     * TODO: refactor this to derive TreeLikelihood!
	 */
    public GraphLikelihood(
            GraphModel graphModel,
            PartitionModel partitionModel,
            TipPartialsModel tipPartialsModel,
            boolean useAmbiguities,
            boolean allowMissingTaxa,
            boolean forceRescaling) 
    {
        super(GRAPH_LIKELIHOOD, createConcatenatedSiteList(partitionModel), graphModel);
        createConcatenatedSiteMap(partitionModel, siteMap);

        try{

        this.concatSiteList = (SiteList)this.patternList;
        this.partitionModel = partitionModel;

        // BEGIN TreeLikelihood LIKELIHOODCORE INIT CODE
        this.likelihoodCore = new GeneralGraphLikelihoodCore(stateCount);
        String coreName = "Java general";
        // END TreeLikelihood LIKELIHOODCORE INIT CODE

        SiteModel siteModel = partitionModel.getPartition(0).getSiteModel();
        this.integrateAcrossCategories = siteModel.integrateAcrossCategories();
        this.categoryCount = siteModel.getCategoryCount();

    	updatePattern = new boolean[patternList.getPatternCount()];
    	updateNode = new boolean[graphModel.getNodeCount()];
    	updatePartition = new boolean[partitionModel.getPartitionCount()];
    	updateNodePartition = new boolean[graphModel.getNodeCount()][partitionModel.getPartitionCount()];

    	
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

	    // mark everything dirty initially
    	updateAllNodes();
	}
    
    protected static SiteList createConcatenatedSiteList(PartitionModel pm){
    	// assume all SiteList instances have the same taxon list.
    	SimpleSiteList ssl = new SimpleSiteList(pm.getSiteList(0).getDataType(),pm.getSiteList(0));
    	for(int i=0; i<pm.getSiteListCount(); i++){
    		SiteList sl = pm.getSiteList(i);
    		for(int j=0; j<sl.getSiteCount(); j++){
        		ssl.addPattern(sl.getSitePattern(j));
    		}
    	}
    	return ssl;
    }

    protected void createConcatenatedSiteMap(PartitionModel pm, HashMap<SiteList, Integer> siteListMap){
    	int sum = 0;
    	for(int i=0; i<pm.getSiteListCount(); i++){
    		SiteList sl = pm.getSiteList(i);
    		siteListMap.put(sl, sum);
    		sum += sl.getSiteCount();
    	}
    }

    
    
    protected void growNodeStorage(){
    	likelihoodCore.growNodeStorage(updateNode.length);	// double the number
    	nodeCount *= 2;

    	// FIXME: need to update all the data structures here
    	boolean[] tmp = new boolean[updateNode.length*2];
    	System.arraycopy(updateNode, 0, tmp, 0, updateNode.length);
    	updateNode = tmp;
    	boolean[][] tmp2 = new boolean[nodeCount][updateNodePartition[0].length];
    	for(int i=0; i<updateNodePartition.length; i++){
    		System.arraycopy(updateNodePartition[i], 0, tmp2[i], 0, tmp2[i].length);
    	}
    	updateNodePartition = tmp2;
    }

    /**
     * Grow the number of partitions the likelihood calc can handle
     * Call this when number of partitions in PartitionModel exceeds current storage capacity
     */
    protected void growPartitionStorage(){
    	int newPartCount = updatePartition.length*2;
    	boolean[] tmp = new boolean[newPartCount];
    	System.arraycopy(updatePartition, 0, tmp, 0, updatePartition.length);
    	updatePartition = tmp;
    	boolean[][] tmp2 = new boolean[nodeCount][newPartCount];
    	for(int i=0; i<nodeCount; i++)
        	System.arraycopy(updateNodePartition[i], 0, tmp2[i], 0, updateNodePartition[i].length);
    	updateNodePartition = tmp2;
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
    	while(treeModel.getNodeCount() >= nodeCount)
    	{
    		growNodeStorage(); // need more node storage!
    	}
    	while(partitionModel.getPartitionCount() > updatePartition.length){
    	    growPartitionStorage();
    	}
    	if(model instanceof PartitionModel){
    		// the partitions changed, do something here
    		// just mark all affected site patterns as dirty
    		if(object instanceof PartitionChangedEvent){
	    		PartitionChangedEvent pce = (PartitionChangedEvent)object;
	    		updatePartition[pce.getPartition().getNumber()] = true;
	    		if(pce.hasNewSection())
	    		{
	    			// mark only the new section as dirty
	    			for(int i=pce.getNewSectionLeft(); i<pce.getNewSectionRight(); i++)
	    			{
	    				updatePattern[i]=true;
	    			}
	    		}
    		}else{
    			// a whole partition is dirty, possibly due to a siteModel or branchRateModel change
    			updatePartition[index] = true;
    		}
    	}else if(model instanceof GraphModel){    		
    		// otherwise the TreeLikelihood can handle it
    		if(object instanceof TreeChangedEvent){
	    		TreeChangedEvent tce = (TreeChangedEvent)object;
	    		updateNode[tce.getNode().getNumber()]=true;
    		}
//    		else if(object instanceof Parameter.Default){
//    			updateNode[((Parameter.Default)object).];
//    		}
    	}

    	super.handleModelChangedEvent(model, object, index);
    }

    protected void storeState() {
        likelihoodCore.storeState();
        super.storeState();
    }
    protected void restoreState() {
        likelihoodCore.restoreState();
        super.restoreState();
    }

    static int ldumpcount =0;
    /**
     * Calculate the log likelihood of the current state.
     * Ripped nearly verbatim from treeLikelihood
     * TODO: refactor this to derive TreeLikelihood!
     * @return the log likelihood.
     */
    protected double calculateLogLikelihood() {

        if (patternLogLikelihoods == null) {
            patternLogLikelihoods = new double[patternCount];
        }

        final NodeRef root = treeModel.getRoot();
        traverse((GraphModel)treeModel, root, partitionModel);

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
            traverse((GraphModel)treeModel, root, partitionModel);

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

        // dump root ll for debug
        ldumpcount++;
        System.err.println("root ll count" + ldumpcount);        

        return logL;
    }
    
    /**
     * partition-aware traversal
     */
    protected void traverse(GraphModel gm, NodeRef node, PartitionModel partitionModel) 
    {
    	if(node==null)	return;	// nothing to see here

        int nodeNum = node.getNumber();

        GraphModel.Node child1 = (GraphModel.Node)gm.getChild(node, 0);
        GraphModel.Node child2 = (GraphModel.Node)gm.getChild(node, 1);

        // If the node has children update the partial likelihoods.
        traverse(gm, child1, partitionModel);
        traverse(gm, child2, partitionModel);

        for(int pI=0; pI < partitionModel.getPartitionCount(); pI++){
        	Partition p = partitionModel.getPartition(pI);

        	// which parent has this partition?
            NodeRef parent = null;
            GraphModel.Node gnode = (GraphModel.Node)node;
            if(gnode.hasObject(0, p))
            	parent = gm.getParent(node, 0);
        	else if(gnode.hasObject(1, p))
        		parent = gm.getParent(node, 1);

            if(parent==null&&gm.getParent(node)!=null){
            	// partition does not pass through this node.
            	continue;
            }

            // does this one need to be updated?
            boolean update = false;
            if (updateNode[nodeNum] || updateNodePartition[nodeNum][p.getNumber()] || updatePartition[p.getNumber()]) {
                update = true;
            }
                    
            // If the node has children update the partial likelihoods.
            if (gm.getChildCount(node)>0 && update)
            {
            	int c1p = gm.getParent(child1, 0)==node ? 0 : 1;
            	int c2p = (child2!=null && gm.getParent(child2, 0)==node) ? 0 : 1;
                boolean has1 = gm.hasObjectOnEdge(child1, c1p, p);
                boolean has2 = child2 != null ? gm.hasObjectOnEdge(child2, c2p, p) : false;
            	// handle reassortment loop corner case
                if(child1==child2) {has1 = false; has2 = false; }

                // determine the left and right bounds of this partition
                int l = remapSite(p.getSiteList(), p.getLeftSite());
                int r = remapSite(p.getSiteList(), p.getRightSite());

                // If we have two children with the partition, and
                // either child node was updated then update this node too
                if (has1&&has2) 
                {
                    // First update the transition probability matrix(ices) for child branches
                	// and get the node number where the partition bifurcates
                    final int peelChild1 = setNodeMatrix(gm, p, child1);
                    final int peelChild2 = setNodeMatrix(gm, p, child2);

                    likelihoodCore.setNodePartialsForUpdate(nodeNum);

                    if (integrateAcrossCategories) {
                        likelihoodCore.calculatePartials(peelChild1, peelChild2, nodeNum, l, r);
                    } else {
                        likelihoodCore.calculatePartials(peelChild1, peelChild2, nodeNum, siteCategories, l, r);
                    }
                    if(l==0)
                    	likelihoodCore.debugPrintPartials(peelChild1, peelChild2, nodeNum, 0);

                    if (COUNT_TOTAL_OPERATIONS) {
                        totalOperationCount ++;
                    }
                }
                if (parent == null) {
                    // No parent this is the root of the tree -
                    // calculate the pattern likelihoods
                	// do this at whatever node the partition coalesed at
                	GraphModel.Node partitionRootNode = getPartitionRoot(p, gm, gnode);
                    FrequencyModel frequencyModel = p.getSiteModel().getFrequencyModel();
                    double[] frequencies = frequencyModel.getFrequencies();
                    double[] partials = getRootPartials(p, partitionRootNode);
                    likelihoodCore.calculateLogLikelihoods(partials, frequencies, patternLogLikelihoods, l, r);
                }
            }
            if(update&&parent!=null)
            	updateNodePartition[parent.getNumber()][p.getNumber()]=true;
        }
    }

    protected GraphModel.Node getPartitionRoot(Partition p, GraphModel gm, GraphModel.Node node){
        GraphModel.Node child1 = (GraphModel.Node)gm.getChild(node, 0);
        GraphModel.Node child2 = (GraphModel.Node)gm.getChild(node, 1);
    	int c1p = gm.getParent(child1, 0)==node ? 0 : 1;
    	int c2p = (child2!=null && gm.getParent(child2, 0)==node) ? 0 : 1;
        boolean has1 = gm.hasObjectOnEdge(child1, c1p, p);
        boolean has2 = child2 != null ? gm.hasObjectOnEdge(child2, c2p, p) : false;
        if(has1&&has2)	return node;
        if(has1&&!has2)	return getPartitionRoot(p, gm, child1);
        if(!has1&&has2)	return getPartitionRoot(p, gm, child2);
        throw new RuntimeException("Error, partition does not trace to root");
    }

    /**
     * Set the child node's matrix immediately prior to integrating partials
     * @param g
     * @param p
     * @param n
     */
    protected int setNodeMatrix(GraphModel g, Partition p, NodeRef n)
    {
        int[] peelChild = new int[1];	// the descendant where the partition bifurcates
        double branchTime = getRateTime(g, p, n, peelChild);
        System.out.print("");
        int nodeNum = n.getNumber();

        likelihoodCore.setNodeMatrixForUpdate(peelChild[0]);

        SiteModel siteModel = p.getSiteModel();
        for (int i = 0; i < categoryCount; i++) {
            double branchLength = siteModel.getRateForCategory(i) * branchTime;
            siteModel.getSubstitutionModel().getTransitionProbabilities(branchLength, probabilities);
            likelihoodCore.setNodeMatrix(peelChild[0], i, probabilities);
        }
        return peelChild[0];
    }
    
    /*
     * computes the total product of rate and time for
     * on the branch at or below n which contains partition p
     * Handles the situation where a partition may be defined at only one
     * child node, so the total rate*time to integrate over includes the
     * next descendant also
     * @param peelChild	an array of size 1, the value will be set to the node id of the
     * 					child where the partition bifurcates
     */
    protected double getRateTime(GraphModel g, Partition p, NodeRef n, int[] peelChild)
    {
    	peelChild[0] = n.getNumber();
        BranchRateModel branchRateModel = p.getBranchRateModel();
        NodeRef parent = g.getParent(n);
    	// find nearest node at or below n where p is
    	// assigned to both children
    	double branchRate = branchRateModel.getBranchRate(g, n);
        // Get the operational time of the branch
        double branchRateTime = branchRate * (g.getNodeHeight(parent) - g.getNodeHeight(n));

        if (branchRateTime < 0.0) {
            throw new RuntimeException("Negative branch length: " + branchRateTime);
        }
        GraphModel.Node c1 = (GraphModel.Node)g.getChild(n, 0);
        GraphModel.Node c2 = (GraphModel.Node)g.getChild(n, 1);
        if(g.getChildCount(n)==1){
        	int c1p = g.getParent(c1, 0)==n ? 0 : 1;
        	if(g.hasObjectOnEdge(c1, c1p, p))
    			branchRateTime += getRateTime(g, p, c1, peelChild);
    		else
    			System.err.println("Partition dead-end");
        }else if(g.getChildCount(n)==2){
        	int c1p = g.getParent(c1, 0)==n ? 0 : 1;
        	int c2p = g.getParent(c2, 0)==n ? 0 : 1;
        	boolean has1 = g.hasObjectOnEdge(c1, c1p, p);
        	boolean has2 = g.hasObjectOnEdge(c2, c2p, p);
            if(c1==c2&&!has1) has1 = g.hasObjectOnEdge(c1, 1, p);
            if(c1==c2&&!has2) has2 = g.hasObjectOnEdge(c1, 1, p);

            if(has1&&!has2)
    			branchRateTime += getRateTime(g, p, g.getChild(n, 0), peelChild);
        	if(!has1&&has2)
    			branchRateTime += getRateTime(g, p, g.getChild(n, 1), peelChild);
        	if(has1&&has2&&c1==c2)	// pick one to follow arbitrarily
    			branchRateTime += getRateTime(g, p, g.getChild(n, 0), peelChild);
        	if(!has1&&!has2)
    			System.err.println("Partition dead-end");
        }
    	return branchRateTime;
    }

    public final double[] getRootPartials(Partition p, NodeRef coalesceNode) {
        if (rootPartials == null) {
            rootPartials = new double[patternCount * stateCount];
        }

        int nodeNum = coalesceNode.getNumber();
        int l = remapSite(p.getSiteList(), p.getLeftSite());
        int r = remapSite(p.getSiteList(), p.getRightSite());
        if (integrateAcrossCategories) {
            // moved this call to here, because non-integrating siteModels don't need to support it - AD
            SiteModel siteModel = p.getSiteModel();
            double[] proportions = siteModel.getCategoryProportions();
            likelihoodCore.integratePartials(nodeNum, proportions, rootPartials, l, r);
        } else {
            likelihoodCore.getPartials(nodeNum, rootPartials, l, r);
        }

        return rootPartials;
    }
    
    /*
     * Converts a partition-local site index to a concatenated site index
     */
    protected int remapSite(SiteList sl, int site){
    	Integer left = siteMap.get(sl);
    	return left.intValue() + site;
    }

    /**
     * the root partial likelihoods (a temporary array that is used
     * to fetch the partials - it should not be examined directly -
     * use getRootPartials() instead).
     */
    private double[] rootPartials = null;
    

    protected PartitionModel partitionModel;
    protected GeneralGraphLikelihoodCore likelihoodCore;
    protected int categoryCount;
    protected final boolean integrateAcrossCategories;
    protected double[] probabilities;
    protected double[] tipPartials;
    protected double[] patternLogLikelihoods = null;
    protected int[] siteCategories = null;
    protected boolean[] updatePartition = null;
    protected boolean[][] updateNodePartition = null;
    
    protected HashMap<SiteList, Integer> siteMap = new HashMap<SiteList, Integer>();
}

/*
 * TreeLikelihood.java
 *
 * (c) 2002-2005 BEAST Development Core Team
 *
 * This package may be distributed under the
 * Lesser Gnu Public Licence (LGPL)
 */

package dr.evomodel.arg.likelihood;

import dr.evolution.alignment.PatternList;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.TaxonList;
import dr.evomodel.arg.ARGModel;
import dr.evomodel.arg.ARGTree;
import dr.evomodel.arg.operators.ARGPartitioningOperator;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.DefaultBranchRateModel;
import dr.evomodel.sitemodel.SiteModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.treelikelihood.*;
import dr.inference.model.Likelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.logging.Logger;
import java.util.*;

/**
 * ARGLikelihood - implements a Likelihood Function for sequences on an ancestral recombination graph.
 *
 * @author Marc Suchard
 * @version $Id: ARGLikelihood.java,v 1.3 2006/10/23 04:13:41 msuchard Exp $
 */

public class ARGLikelihood extends AbstractARGLikelihood {

	public static final String ARG_LIKELIHOOD = "argTreeLikelihood";
	public static final String USE_AMBIGUITIES = "useAmbiguities";
	public static final String STORE_PARTIALS = "storePartials";
	public static final String USE_SCALING = "useScaling";

	/**
	 * Constructor.
	 */
	public ARGLikelihood(PatternList patternList,
	                     ARGModel treeModel,
	                     SiteModel siteModel,
	                     BranchRateModel branchRateModel,
	                     boolean useAmbiguities,
	                     boolean storePartials,
	                     boolean useScaling) {

		super(ARG_LIKELIHOOD, patternList, treeModel);

		partition = treeModel.addLikelihoodCalculator(this);
		
		this.storePartials = storePartials;
		this.useAmbiguities = useAmbiguities;

		try {
			this.siteModel = siteModel;
			addModel(siteModel);

			this.frequencyModel = siteModel.getFrequencyModel();
			addModel(frequencyModel);

			integrateAcrossCategories = siteModel.integrateAcrossCategories();

			this.categoryCount = siteModel.getCategoryCount();

			if (integrateAcrossCategories) {
				if (patternList.getDataType() instanceof dr.evolution.datatype.Nucleotides) {

					if (NativeNucleotideLikelihoodCore.isAvailable()) {

						Logger.getLogger("dr.evomodel").info("TreeLikelihood using native nucleotide likelihood core");
						likelihoodCore = new NativeNucleotideLikelihoodCore();
					} else {

						Logger.getLogger("dr.evomodel").info("TreeLikelihood using Java nucleotide likelihood core");
						likelihoodCore = new NucleotideLikelihoodCore();
					}

				} else if (patternList.getDataType() instanceof dr.evolution.datatype.AminoAcids) {

					if (NativeAminoAcidLikelihoodCore.isAvailable()) {
						Logger.getLogger("dr.evomodel").info("TreeLikelihood using native amino acid likelihood core");
						likelihoodCore = new NativeAminoAcidLikelihoodCore();
					} else {
						Logger.getLogger("dr.evomodel").info("TreeLikelihood using java likelihood core");
						likelihoodCore = new AminoAcidLikelihoodCore();
					}
				} else if (patternList.getDataType() instanceof dr.evolution.datatype.Codons) {
					Logger.getLogger("dr.evomodel").info("TreeLikelihood using Java codon likelihood core");
//					likelihoodCore = new CodonLikelihoodCore(patternList.getStateCount());
					this.useAmbiguities = true;
					throw new RuntimeException("Still need to merge codon likelihood core");
				} else {
					if (patternList.getDataType() instanceof dr.evolution.datatype.HiddenNucleotides &&
							NativeCovarionLikelihoodCore.isAvailable()) {
						Logger.getLogger("dr.evomodel").info("TreeLikelihood using native covarion likelihood core");
						likelihoodCore = new NativeCovarionLikelihoodCore();
					} else {
						Logger.getLogger("dr.evomodel").info("TreeLikelihood using Java general likelihood core");
						likelihoodCore = new GeneralLikelihoodCore(patternList.getStateCount());
					}

				}
			} else {

				Logger.getLogger("dr.evomodel").info("TreeLikelihood using Java general likelihood core");
				likelihoodCore = new GeneralLikelihoodCore(patternList.getStateCount());
			}
			Logger.getLogger("dr.evomodel").info("  " + (useAmbiguities ? "Using" : "Ignoring") + " ambiguities in tree likelihood.");
			Logger.getLogger("dr.evomodel").info("  Partial likelihood scaling " + (useScaling ? "on." : "off."));

			if (branchRateModel != null) {
				this.branchRateModel = branchRateModel;
				Logger.getLogger("dr.evomodel").info("Branch rate model used: " + branchRateModel.getModelName());
			} else {
				this.branchRateModel = new DefaultBranchRateModel();
			}
			addModel(this.branchRateModel);

			probabilities = new double[stateCount * stateCount];

//			likelihoodCore.initialize(nodeCount, patternCount, categoryCount, integrateAcrossCategories, useScaling);
			likelihoodCore.initialize(nodeCount, patternCount, categoryCount, integrateAcrossCategories);

			int extNodeCount = treeModel.getExternalNodeCount();
			int intNodeCount = treeModel.getInternalNodeCount();

			for (int i = 0; i < extNodeCount; i++) {
				// Find the id of tip i in the patternList
				String id = treeModel.getTaxonId(i);
				int index = patternList.getTaxonIndex(id);

//                System.err.println("id = "+id+"  index = "+index);

				if (index == -1) {
					throw new TaxonList.MissingTaxonException("Taxon, " + id + ", in tree, " + treeModel.getId() +
							", is not found in patternList, " + patternList.getId());
				}

				if (useAmbiguities) {
					setPartials(likelihoodCore, patternList, categoryCount, index, i);
				} else {
					setStates(likelihoodCore, patternList, index, i);
				}
			}

//            System.exit(-1);
			for (int i = 0; i < intNodeCount; i++) {
				likelihoodCore.createNodePartials(extNodeCount + i);
			}
		} catch (TaxonList.MissingTaxonException mte) {
			throw new RuntimeException(mte.toString());
		}

	}

	// **************************************************************
	// ModelListener IMPLEMENTATION
	// **************************************************************


    private static final boolean NO_CACHING = false;

    /**
     * Handles model changed events from the submodels.
     */
    protected void handleModelChangedEvent(Model model, Object object, int index) {

        if (NO_CACHING) {
            reconstructTree = true;
            updateAllNodes();
        }

        if (model == treeModel) {
            if (object instanceof ARGModel.TreeChangedEvent) {
                ARGModel.TreeChangedEvent event = (ARGModel.TreeChangedEvent) object;
                if (event.isSizeChanged() ) {
                    updateAllNodes(); // TODO Update only affected portion of tree
                    reconstructTree = true;
                } else if (event.isNodeChanged()) {
                    // If a node event occurs the node and its two child nodes
                    // are flagged for updating (this will result in everything
                    // above being updated as well. Node events occur when a node
                    // is added to a branch, removed from a branch or its height or
                    // rate changes.

                    NodeRef treeNode = mapARGNodesToTreeNodes.get(event.getNode());
                    if ( treeNode != null ) {                        
                        if (event.isHeightChanged() || event.isRateChanged()) {
                            updateNodeAndChildren(treeNode);
                        } else {
                            reconstructTree = true;
//                            updateNodeAndChildren(treeNode); // TODO This doesn't work with sizeChange; why???
                            updateAllNodes();
                        }
                    } 
                } else if (event.isTreeChanged()) {
                    // Full tree events result in a complete updating of the tree likelihood
                    // These include adding and removing nodes
                    // TODO ARG rearrangements still call this; they should not      
                    reconstructTree = true;
                    updateAllNodes();
                } else {
                    // Other event types are ignored (probably trait changes).
                    throw new RuntimeException("Another tree event has occured (possibly a trait change).");
                }
            } else if (object instanceof ARGPartitioningOperator.PartitionChangedEvent) {
	      final boolean[] updatePartition = ((ARGPartitioningOperator.PartitionChangedEvent) object).getUpdatedPartitions();
                if (updatePartition[partition]) {
                    reconstructTree = true;
                    updateAllNodes(); // TODO Probably does not affect entire tree; fix
                }
            } else if (object instanceof Parameter) {
                // ignore, most of these are handled in isNodeChanged()
            } else
                throw new RuntimeException("Unexpected ARGModel update "+object.getClass());

        } else if (model == branchRateModel) {
            // TODO Only update affected branches
            updateAllNodes();
        } else if (model == frequencyModel) {
            updateAllNodes();
        } else if (model instanceof SiteModel) {
            updateAllNodes();
        } else {
            throw new RuntimeException("Unknown componentChangedEvent");
        }

        super.handleModelChangedEvent(model, object, index);
    }

	// **************************************************************
	// Model IMPLEMENTATION
	// **************************************************************

	/**
	 * Stores the additional state other than model components
	 */
          protected void storeState() {

              if (storePartials) {
                  likelihoodCore.storeState();
              }
              super.storeState();
          }

	/**
	 * Restore the additional stored state
	 */
          protected void restoreState() {

              if (storePartials) {
                  likelihoodCore.restoreState();
              } else {
                  updateAllNodes();
              }
              reconstructTree = true; // currently the tree is not cached, because the ARG that generates it is cached
              super.restoreState();
          }

        private int getUnusedInt(Map<NodeRef,Integer> inMap) {
            Collection<Integer> intSet = inMap.values();
            int i = tree.getExternalNodeCount();
            while( intSet.contains(i) )
                i++;
            return i;
        }


    private Set<NodeRef> unsetNodes = null;

        private void reconstructTree() {
            
            oldTree = tree;
            oldMapARGNodesToInts = mapARGNodesToInts;

            tree = new ARGTree(treeModel, partition);            
            reconstructTree = false;
            mapARGNodesToInts = new HashMap<NodeRef,Integer>(tree.getInternalNodeCount());
            mapARGNodesToTreeNodes = tree.getMapping();

            if (oldTree == null) {
                 // First initialization
                for(int i=0; i<tree.getInternalNodeCount(); i++) {
                    NodeRef node = tree.getInternalNode(i);
                    mapARGNodesToInts.put(treeModel.getMirrorNode(node),node.getNumber());
                }
            } else {

                // Need to renumber
                if (unsetNodes == null)
                    unsetNodes = new HashSet<NodeRef>();
                else
                    unsetNodes.clear();

                // Copy over numbers for nodes that still exist in tree
                for (int i = 0; i < tree.getInternalNodeCount(); i++) {
                    NodeRef newNode = tree.getInternalNode(i);
                    NodeRef argNode = treeModel.getMirrorNode(newNode);

                    if (oldMapARGNodesToInts.containsKey(argNode)) { // was in old tree

                        int oldNumber = oldMapARGNodesToInts.get(argNode);
                        treeModel.setNodeNumber(newNode,oldNumber);
                        mapARGNodesToInts.put(argNode,oldNumber);
                    } else  // was not in old tree
                        unsetNodes.add(newNode);
                }

                // Set unused numbers for nodes that are new and mark for update
                for (NodeRef node : unsetNodes) {
                    int newNumber = getUnusedInt(mapARGNodesToInts);
                    treeModel.setNodeNumber(node,newNumber);
                    mapARGNodesToInts.put(node,newNumber);
                    updateNode[newNumber] = true;
                }
            }
        }

	// **************************************************************
	// Likelihood IMPLEMENTATION
	// **************************************************************

	/**
	 * Calculate the log likelihood of the current state.
	 *
	 * @return the log likelihood.
	 */
	protected double calculateLogLikelihood() {

              if (reconstructTree) {
                  reconstructTree();
              }
              
		NodeRef root = tree.getRoot();

		if (rootPartials == null) {
			rootPartials = new double[patternCount * stateCount];
		}

		if (patternLogLikelihoods == null) {
			patternLogLikelihoods = new double[patternCount];
		}

		if (!integrateAcrossCategories) {
			if (siteCategories == null) {
				siteCategories = new int[patternCount];
			}
			for (int i = 0; i < patternCount; i++) {
				siteCategories[i] = siteModel.getCategoryOfSite(i);
			}
		}

        try {
		    traverse(tree, root);
        } catch (NegativeBranchLengthException e) {
            System.err.println("Negative branch length found, trying to return 0 likelihood");
            return Double.NEGATIVE_INFINITY;
        }

		//********************************************************************
		// after traverse all nodes and patterns have been updated --
		//so change flags to reflect this.
		for (int i = 0; i < nodeCount; i++) {
			updateNode[i] = false;
		}
		//********************************************************************

		double logL = 0.0;

		for (int i = 0; i < patternCount; i++) {
//        	System.err.printf("Pattern %2d:  %5.4f  %5.4f\n",i,patternLogLikelihoods[i],patternWeights[i]);
			logL += patternLogLikelihoods[i] * patternWeights[i];
		}

		return logL;
	}

    class NegativeBranchLengthException extends Exception {

    }

	/**
	 * Traverse the tree calculating partial likelihoods.
	 *
	 * @return whether the partials for this node were recalculated.
	 */
	private boolean traverse(Tree tree, NodeRef node) throws NegativeBranchLengthException {

		boolean update = false;

		int nodeNum = node.getNumber();
//        System.err.println(nodeNum);

		NodeRef parent = tree.getParent(node);
		
		

		// First update the transition probability matrix(ices) for this branch
		if (parent != null && updateNode[nodeNum]) {


			double branchRate = branchRateModel.getBranchRate(tree, node);
						
			// Get the operational time of the branch
			double branchTime = branchRate * (tree.getNodeHeight(parent) - tree.getNodeHeight(node));
			if (branchTime < 0.0) {
                if (!DEBUG) {
				    throw new RuntimeException("Negative branch length: " + branchTime);
                } else{
                    throw new NegativeBranchLengthException();
                }
			}

			for (int i = 0; i < categoryCount; i++) {
                double branchLength = siteModel.getRateForCategory(i) * branchTime;
                siteModel.getSubstitutionModel().getTransitionProbabilities(branchLength, probabilities);
				likelihoodCore.setNodeMatrix(nodeNum, i, probabilities);
			}

			update = true;
		}

		// If the node is internal, update the partial likelihoods.
		if (!tree.isExternal(node)) {

			int nodeCount = tree.getChildCount(node);
			if (nodeCount != 2)
				throw new RuntimeException("binary trees only!");

			// Traverse down the two child nodes
			NodeRef child1 = tree.getChild(node, 0);
			boolean update1 = traverse(tree, child1);

			NodeRef child2 = tree.getChild(node, 1);
			boolean update2 = traverse(tree, child2);

			// If either child node was updated then update this node too
			if (update1 || update2) {

				int childNum1 = child1.getNumber();
				int childNum2 = child2.getNumber();

				if (integrateAcrossCategories) {
					likelihoodCore.calculatePartials(childNum1, childNum2, nodeNum);
				} else {
					likelihoodCore.calculatePartials(childNum1, childNum2, nodeNum,
							siteCategories);
				}

				if (parent == null) {
					// No parent this is the root of the tree -
					// calculate the pattern likelihoods
					double[] frequencies = frequencyModel.getFrequencies();

					if (integrateAcrossCategories) {

						// moved this call to here, because non-integrating siteModels don't need to support it - AD
						double[] proportions = siteModel.getCategoryProportions();
						likelihoodCore.integratePartials(nodeNum, proportions, rootPartials);
					} else {
						likelihoodCore.getPartials(nodeNum, rootPartials);
					}

					likelihoodCore.calculateLogLikelihoods(rootPartials, frequencies, patternLogLikelihoods);
				}

				update = true;
			}
		}

		return update;

	}

	/**
	 * The XML parser
	 */
	public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

		public String getParserName() {
			return ARG_LIKELIHOOD;
		}

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {

			boolean useAmbiguities = false;
			boolean storePartials = true;
			boolean useScaling = false;
			if (xo.hasAttribute(USE_AMBIGUITIES)) {
				useAmbiguities = xo.getBooleanAttribute(USE_AMBIGUITIES);
			}
			if (xo.hasAttribute(STORE_PARTIALS)) {
				storePartials = xo.getBooleanAttribute(STORE_PARTIALS);
			}
			if (xo.hasAttribute(USE_SCALING)) {
				useScaling = xo.getBooleanAttribute(USE_SCALING);
			}
			PatternList patternList = (PatternList) xo.getChild(PatternList.class);
			ARGModel treeModel = (ARGModel) xo.getChild(ARGModel.class);
			SiteModel siteModel = (SiteModel) xo.getChild(SiteModel.class);

			BranchRateModel branchRateModel = (BranchRateModel) xo.getChild(BranchRateModel.class);

			return new ARGLikelihood(patternList, treeModel, siteModel, branchRateModel, useAmbiguities, storePartials, useScaling);
		}

		//************************************************************************
		// AbstractXMLObjectParser implementation
		//************************************************************************

		public String getParserDescription() {
			return "This element represents the likelihood of a patternlist on a tree given the site model.";
		}

		public Class getReturnType() {
			return Likelihood.class;
		}

		public XMLSyntaxRule[] getSyntaxRules() {
			return rules;
		}

		private final XMLSyntaxRule[] rules = {
				AttributeRule.newBooleanRule(USE_AMBIGUITIES, true),
				AttributeRule.newBooleanRule(STORE_PARTIALS, true),
				AttributeRule.newBooleanRule(USE_SCALING, true),
				new ElementRule(PatternList.class),
				new ElementRule(ARGModel.class),
				new ElementRule(SiteModel.class),
				new ElementRule(BranchRateModel.class, true)
		};
	};

	/**
	 * XML Serializer for parallelization
	 *
	 */

//	public Element toXML() {
//		Element likelihoodElement
//	}

	// **************************************************************
	// INSTANCE VARIABLES
	// **************************************************************

	/**
	 * the frequency model for these sites
	 */
	protected FrequencyModel frequencyModel = null;

	/**
	 * the site model for these sites
	 */
	protected SiteModel siteModel = null;

	/**
	 * the branch rate model
	 */
	protected BranchRateModel branchRateModel = null;

	private boolean storePartials = false;

	private boolean integrateAcrossCategories = false;

	/**
	 * the categories for each site
	 */
	protected int[] siteCategories = null;

	/**
	 * the root partial likelihoods
	 */
	protected double[] rootPartials = null;

	/**
	 * the pattern likelihoods
	 */
	protected double[] patternLogLikelihoods = null;

	/**
	 * the number of rate categories
	 */
	protected int categoryCount;

	/**
	 * an array used to store transition probabilities
	 */
	protected double[] probabilities;

	/**
	 * the LikelihoodCore
	 */
	protected LikelihoodCore likelihoodCore;

	private boolean useAmbiguities;

    private boolean reconstructTree = true;
    private ARGTree tree = null;
    private ARGTree oldTree;

    private Map<NodeRef,Integer> mapARGNodesToInts = null;
    private Map<NodeRef,Integer> oldMapARGNodesToInts;

    private Map<NodeRef,NodeRef> mapARGNodesToTreeNodes = null;

    private static final boolean DEBUG = true;
}
/*
 * TipsTreeLikelihood.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.evomodel.treelikelihood; 

import dr.evolution.alignment.PatternList;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.TaxonList;
import dr.evomodel.sitemodel.SiteModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Likelihood;
import dr.inference.model.Model;
import dr.xml.*;

/**
 * TipTreeLikelihood - implements a Likelihood Function for sequences on a tree.
 * This one has allows a different siteModel at the end of each tip to represent
 * DNA damage or other Tip-Influencing-Processes.
 *
 * @version $Id: TipsTreeLikelihood.java,v 1.6 2006/01/10 16:48:28 rambaut Exp $
 *
 * @author Andrew Rambaut
 */

public class TipsTreeLikelihood extends AbstractTreeLikelihood {

	public static final String TIPS_TREE_LIKELIHOOD = "tipsTreeLikelihood";
	public static final String TIPS = "tips";

    /**
     * Constructor.
     */   	
    public TipsTreeLikelihood(	PatternList patternList,
    								TreeModel treeModel,
    								SiteModel siteModel,
    								SiteModel tipsSiteModel,
                                    boolean useScaling) throws TaxonList.MissingTaxonException
	{
    	
		super(TIPS_TREE_LIKELIHOOD, patternList, treeModel);
		
		try {
			this.siteModel = siteModel;
			addModel(siteModel);
			
			this.tipsSiteModel = tipsSiteModel;
			addModel(tipsSiteModel);
			
			this.frequencyModel = siteModel.getFrequencyModel();
			addModel(frequencyModel);
	    	
	    	if (!siteModel.integrateAcrossCategories()) {
	    		throw new RuntimeException("TipsTreeLikelihood can only use SiteModels that require integration across categories");
	    	}
	    	
	 		this.categoryCount = siteModel.getCategoryCount();
	   		   		
			if (patternList.getDataType() instanceof dr.evolution.datatype.Nucleotides) {

				if (NativeNucleotideLikelihoodCore.isAvailable()) {
				
					System.out.println("TipsTreeLikelihood using native nucleotide likelihood core.");
					likelihoodCore = new NativeNucleotideLikelihoodCore();
					tipsLikelihoodCore = new NativeNucleotideLikelihoodCore();
				} else {
				
					System.out.println("TipsTreeLikelihood using Java nucleotide likelihood core.");
					likelihoodCore = new NucleotideLikelihoodCore();
					tipsLikelihoodCore = new NucleotideLikelihoodCore();
				}
				
			} else if (patternList.getDataType() instanceof dr.evolution.datatype.AminoAcids) {
				System.out.println("TipsTreeLikelihood using Java amino acid likelihood core.");
				likelihoodCore = new AminoAcidLikelihoodCore();
				tipsLikelihoodCore = new AminoAcidLikelihoodCore();
			} else {
				System.out.println("TipsTreeLikelihood using Java general likelihood core.");
				likelihoodCore = new GeneralLikelihoodCore(patternList.getStateCount());
				tipsLikelihoodCore = new GeneralLikelihoodCore(patternList.getStateCount());
			}
			
			probabilities = new double[stateCount * stateCount];

			externalNodeCount = treeModel.getExternalNodeCount();
			int intNodeCount = treeModel.getInternalNodeCount();

			likelihoodCore.initialize(nodeCount, patternCount, categoryCount, true);
			
			// we only need nodes for the tips + 2 extra: 1 for a set of dummy unknown states
			// and one for the temporary destination partials
			tipsLikelihoodCore.initialize(externalNodeCount + 2, patternCount, categoryCount, true);
			
			for (int i = 0; i < externalNodeCount; i++) {
				// Find the id of tip i in the patternList
				String id = treeModel.getTaxonId(i);
				int index = patternList.getTaxonIndex(id);
				
				if (index == -1) {
					throw new TaxonList.MissingTaxonException("Taxon, " + id + ", in tree, " + treeModel.getId() + 
																", is not found in patternList, " + patternList.getId());
				}
				
				// The tipsLikelihoodCore takes the actual data...
				setStates(tipsLikelihoodCore, patternList, index, i);
				
				//... and will set the normal likelihood core tip partials but
				// for the moment we will set them up as normal to allocate the 
				// memory
				setPartials(likelihoodCore, patternList, categoryCount, index, i);
			}
			
			for (int i = 0; i < intNodeCount; i++) {
				likelihoodCore.createNodePartials(externalNodeCount + i);
			}

			// set up the dummy unknown states
			int[] states = new int[patternCount];
			int n = patternList.getDataType().getUnknownState();
			for (int i = 0; i < patternCount; i++) {
				states[i] = n;
			}
			tipsLikelihoodCore.setNodeStates(externalNodeCount, states);

			// and the destination partials
			tipsLikelihoodCore.createNodePartials(externalNodeCount + 1);
			
		} catch (TaxonList.MissingTaxonException mte) {
			throw new RuntimeException(mte.toString());
		}

   	}
    
	// **************************************************************
    // ModelListener IMPLEMENTATION
    // **************************************************************
	
	/**
	 * Handles model changed events from the submodels. 
	 */
	protected void handleModelChangedEvent(Model model, Object object, int index) {

		if (model == treeModel) {
            if (object instanceof TreeModel.TreeChangedEvent) {

                if (((TreeModel.TreeChangedEvent)object).isNodeChanged()) {

                    updateNodeAndChildren(((TreeModel.TreeChangedEvent)object).getNode());

                } else {
                    updateAllNodes();
				
                }
            }
		} else if (model == frequencyModel) {
			
			updateTips = true;
			updateAllNodes();
				
		} else if (model instanceof SiteModel) {
			
			if (model == siteModel) {
				
				updateAllNodes();
					
			} else if (model == tipsSiteModel) {
			
				updateTips = true;
				updateAllNodes();
			}
						
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
		
		likelihoodCore.storeState();
		tipsLikelihoodCore.storeState();
		super.storeState();
		
	}
	
	/**
	 * Restore the additional stored state
	 */
	protected void restoreState() {
	
		likelihoodCore.restoreState();
		tipsLikelihoodCore.restoreState();
		super.restoreState();
			
	}
	
    // **************************************************************
    // Likelihood IMPLEMENTATION
    // **************************************************************

	/**
     * Calculate the log likelihood of the current state.
     * @return the log likelihood.
     */
    protected double calculateLogLikelihood() {
	 
		NodeRef root = treeModel.getRoot();
					
		if (rootPartials == null) {
			rootPartials = new double[patternCount * stateCount];
		}
		
		if (patternLogLikelihoods == null) {
			patternLogLikelihoods = new double[patternCount];
		}
		
		if (updateTips) {
			double[] proportions = tipsSiteModel.getCategoryProportions();
					
			for (int i = 0; i < externalNodeCount; i++) {
		   		for (int j = 0; j < categoryCount; j++) {
					
					tipsSiteModel.getTransitionProbabilitiesForCategory(j, 1.0, probabilities);
					tipsLikelihoodCore.setNodeMatrix(i, j, probabilities);
				}

				tipsLikelihoodCore.calculatePartials(i, externalNodeCount, externalNodeCount + 1);
				
				tipsLikelihoodCore.integratePartials(externalNodeCount + 1, proportions, rootPartials);
				likelihoodCore.setNodePartials(i, rootPartials);
		
			}
			updateTips = false;
		}

		traverse(treeModel, root);
		
		//********************************************************************
		// after traverse all nodes and patterns have been updated -- 
		//so change flags to reflect this.
		for (int i = 0; i < nodeCount; i++) {
			updateNode[i] = false;
		}
		//********************************************************************
		
		double logL = 0.0;
		
		for (int i = 0; i < patternCount; i++) {
			logL += patternLogLikelihoods[i] * patternWeights[i];
		}
		
        if (Double.isNaN(logL)) {
            throw new RuntimeException("Likelihood NaN");
        }

		return logL;
    }

	/**
     * Traverse the tree calculating partial likelihoods.
     * @return whether the partials for this node were recalculated.
     */
	private final boolean traverse(Tree tree, NodeRef node) {
	
		boolean update = false;
		
		int nodeNum = node.getNumber();
				
		NodeRef parent = tree.getParent(node);
		
		// First update the transition probability matrix(ices) for this branch
		if (parent != null && updateNode[nodeNum]) {

			// Get the average rate over this branch			
	   	   	
	   	   	// ***************************************************************
			// Rate at nodes model
	   		//double branchRate = (tree.getNodeRate(node) + tree.getNodeRate(parent)) / 2;
	   		// ***************************************************************
			
   			// ***************************************************************
	   		// Rate at branches model
	   		double branchRate = tree.getNodeRate(node);
	   		// ***************************************************************
	   		
	   		if (branchRate < 0.0) {
	   			throw new RuntimeException("Negative branch rate: " + branchRate);
	   		}
	   						
			// Get the operational time of the branch
	   		double branchTime = branchRate * ( tree.getNodeHeight(parent) - tree.getNodeHeight(node) );
	   		if (branchTime < 0.0) {
	   			throw new RuntimeException("Negative branch length: " + branchTime);
	   		}
	   		
	   		for (int i = 0; i < categoryCount; i++) {
				
				siteModel.getTransitionProbabilitiesForCategory(i, branchTime, probabilities);
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
			
				likelihoodCore.calculatePartials(childNum1, childNum2, nodeNum);
				
				if (parent == null) {
					// No parent this is the root of the tree - 
					// calculate the pattern likelihoods
					double[] frequencies = frequencyModel.getFrequencies();						
					double[] proportions = siteModel.getCategoryProportions();
					
					likelihoodCore.integratePartials(nodeNum, proportions, rootPartials);
					likelihoodCore.calculateLogLikelihoods(rootPartials, frequencies, patternLogLikelihoods);
				}
																		
				update = true;					
			}
		}
				
		return update;

	}
	
    // **************************************************************
    // XMLElement IMPLEMENTATION
    // **************************************************************

	public org.w3c.dom.Element createElement(org.w3c.dom.Document d) {
		throw new RuntimeException("createElement not implemented");
	}


	/** 
	 * The XML parser
	 */
	public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {
				
		public String getParserName() { return TIPS_TREE_LIKELIHOOD; }
	
		public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            boolean useScaling = false;

            PatternList patternList = (PatternList)xo.getChild(PatternList.class);
			TreeModel treeModel = (TreeModel)xo.getChild(TreeModel.class);
			SiteModel siteModel = (SiteModel)xo.getChild(SiteModel.class);
			
			XMLObject xoc = (XMLObject)xo.getChild(TIPS);
			SiteModel tipsSiteModel = (SiteModel)xoc.getChild(SiteModel.class);
			
			TipsTreeLikelihood treeLikelihood = null;
			try {
				treeLikelihood = new TipsTreeLikelihood(patternList, treeModel, siteModel, tipsSiteModel, useScaling);
			} catch (TaxonList.MissingTaxonException e) {
				throw new XMLParseException(e.toString()); 
			}
			
			return treeLikelihood;
		}
		
		//************************************************************************
		// AbstractXMLObjectParser implementation
		//************************************************************************
		
		public String getParserDescription() {
			return "This element represents the likelihood of a patternlist on a tree given the site model.";
		}
		
		public Class getReturnType() { return Likelihood.class; }
		
		public XMLSyntaxRule[] getSyntaxRules() { return rules; }
		
		private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
			new ElementRule(TIPS,
				new XMLSyntaxRule[] {
					new ElementRule(SiteModel.class, "A siteModel that will be applied only to this set of tips")
				}),
//			new ElementRule(TIPS,
//				new XMLSyntaxRule[] {
//					new ElementRule(Taxa.class, "An optional set of taxa which defines which tips to apply a different site model to. If missing this will be applied to all tips"),
//					new ElementRule(SiteModel.class, "A siteModel that will be applied only to this set of tips")
//				}, 0, Integer.MAX_VALUE),
			new ElementRule(PatternList.class),
			new ElementRule(TreeModel.class),
			new ElementRule(SiteModel.class)
		};
	};
	

	// **************************************************************
	// INSTANCE VARIABLES
	// **************************************************************

	/** the frequency model for these sites */
	protected FrequencyModel frequencyModel = null;

	/** the site model for these sites */
	protected SiteModel siteModel = null;
	
	/** the site models for the tips */
	protected SiteModel tipsSiteModel = null;
	
	private int externalNodeCount;
	
	private boolean updateTips = false;
	    			
	/** the root partial likelihoods */
    protected double[] rootPartials = null;

	/** the pattern likelihoods */
    protected double[] patternLogLikelihoods = null;

	/** the number of rate categories */
	protected int categoryCount;
		
	/** an array used to store transition probabilities */
	protected double[] probabilities;
	
   /** the LikelihoodCore */
	protected LikelihoodCore likelihoodCore;
	
   /** the LikelihoodCore for tips */
	protected LikelihoodCore tipsLikelihoodCore;
}
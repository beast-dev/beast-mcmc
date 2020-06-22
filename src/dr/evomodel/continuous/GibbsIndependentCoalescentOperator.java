/*
 * GibbsIndependentCoalescentOperator.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.continuous;

import dr.evolution.tree.Tree;
import dr.evolution.util.Taxa;
import dr.evolution.util.TaxonList;
import dr.evomodel.coalescent.CoalescentSimulator;
import dr.evomodel.coalescent.DemographicModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.operators.GibbsOperator;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.SimpleMCMCOperator;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

/**
 * An independent coalescent sampler, based on the coalescent simulator parser code.
 *
 * @author Guy Baele
 * 
 */
public class GibbsIndependentCoalescentOperator extends SimpleMCMCOperator implements GibbsOperator {

	public static final String OPERATOR_NAME = "GibbsIndependentCoalescentOperator";
	public static final String HEIGHT = "height";

	private XMLObject xo;
	private TreeModel treeModel;
    private DemographicModel demoModel;
	private double rootHeight;

    private final dr.evolution.coalescent.CoalescentSimulator simulator = new dr.evolution.coalescent.CoalescentSimulator();

	public GibbsIndependentCoalescentOperator(XMLObject xo, TreeModel treeModel, DemographicModel model, double rootHeight, double weight) {

		this.xo = xo;
		this.treeModel = treeModel;
		this.demoModel = model;
		this.rootHeight = rootHeight;
		setWeight(weight);
		
	}

	@Override
	public void setPathParameter(double beta) {
		//do nothing
	}
	
	public String getPerformanceSuggestion() {
		return "";
	}

	public String getOperatorName() {
		return "GibbsIndependentCoalescent(" + treeModel.getModelName() + ")";
	}
	
	public int getStepCount() {
        return 1;
    }

    /**
	 * change the parameter and return the hastings ratio.
     */
	public double doOperation() {

		CoalescentSimulator simulator = new CoalescentSimulator();

		//store the generated tree
		Tree tree;

		List<TaxonList> taxonLists = new ArrayList<TaxonList>();
		List<Tree> subtrees = new ArrayList<Tree>();

		// should have one child that is node
		for (int i = 0; i < xo.getChildCount(); i++) {
			final Object child = xo.getChild(i);

			// AER - swapped the order of these round because Trees are TaxonLists...
			if (child instanceof TreeModel) {
				//do nothing
			} else if (child instanceof Tree) {
				//don't do this as all subtrees will remain identical
				//subtrees.add((Tree) child);

				Tree currentTree = (Tree)child;
				Tree[] newTrees = new Tree[currentTree.getTaxonCount()];
				for (int j = 0; j < currentTree.getTaxonCount(); j++) {
					Taxa tip = new Taxa();
					tip.addTaxon(currentTree.getTaxon(j));
					newTrees[j] = simulator.simulateTree(tip, demoModel);
				}
				Tree newTree = simulator.simulateTree(newTrees, demoModel, rootHeight, newTrees.length != 1);
				subtrees.add(newTree);
			} else if (child instanceof TaxonList) {
				taxonLists.add((TaxonList) child);
			}
		}

		if (taxonLists.size() == 0) {
			if (subtrees.size() == 1) {
				//System.out.println(" *** 1 *** ");
				tree = subtrees.get(0);
				//this would be the normal way to do it
				treeModel.beginTreeEdit();
				//now it's allowed to adjust the tree structure
				treeModel.adoptTreeStructure(tree);
				//endTreeEdit() would then fire the events
				treeModel.endTreeEdit();
				return 0;
			}
			throw new RuntimeException("Expected at least one taxonList or two subtrees in "
					+ getOperatorName() + " XML specification.");
		}

		Taxa remainingTaxa = new Taxa();
		for (int i = 0; i < taxonLists.size(); i++) {
			remainingTaxa.addTaxa(taxonLists.get(i));
		}

		for (int i = 0; i < subtrees.size(); i++) {
			remainingTaxa.removeTaxa(subtrees.get(i));
		}

		Tree[] trees = new Tree[subtrees.size() + remainingTaxa.getTaxonCount()];
		// add the preset trees
		for (int i = 0; i < subtrees.size(); i++) {
			trees[i] = subtrees.get(i);
		}

		// add all the remaining taxa in as single tip trees...
		for (int i = 0; i < remainingTaxa.getTaxonCount(); i++) {
			Taxa tip = new Taxa();
			tip.addTaxon(remainingTaxa.getTaxon(i));
			trees[i + subtrees.size()] = simulator.simulateTree(tip, demoModel);
		}

		tree = simulator.simulateTree(trees, demoModel, rootHeight, trees.length != 1);

		//this would be the normal way to do it
		treeModel.beginTreeEdit();
		//now it's allowed to adjust the tree structure
		treeModel.adoptTreeStructure(tree);
		//endTreeEdit() would then fire the events
		treeModel.endTreeEdit();
		
        return 0;
	}
	
	public static dr.xml.XMLObjectParser PARSER = new dr.xml.AbstractXMLObjectParser() {
		
		public String getParserName() {
            return OPERATOR_NAME;
        }

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {

			double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
			TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
			DemographicModel demoModel = (DemographicModel) xo.getChild(DemographicModel.class);
			double height = xo.getAttribute(HEIGHT, Double.NaN);

			//return new GibbsIndependentCoalescentOperator(trees, treeModel, demoModel, height, trees.length != 1, weight);
			return new GibbsIndependentCoalescentOperator(xo, treeModel, demoModel, height, weight);

		}
		
		//************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

		public XMLSyntaxRule[] getSyntaxRules() {
			return rules;
		}
		
		private final XMLSyntaxRule[] rules = {
				AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
				AttributeRule.newDoubleRule(HEIGHT, true, ""),
                new ElementRule(TreeModel.class),
                new ElementRule(DemographicModel.class),
				new ElementRule(Tree.class, 0, Integer.MAX_VALUE),
				new ElementRule(TaxonList.class, 0, Integer.MAX_VALUE)
        };

		public String getParserDescription() {
			return "This element returns an independence coalescent sampler, disguised as a Gibbs operator, from a demographic model.";
		}

		public Class getReturnType() {
			return MCMCOperator.class;
		}
		
	};

}

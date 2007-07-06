/*
 * ParsimonyStateStatistic.java
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

package dr.evomodel.tree;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxa;
import dr.evolution.util.TaxonList;
import dr.inference.model.Statistic;
import dr.xml.*;

import java.util.Set;


/**
 * A statistic that reconstructs the parsimony state at a mrca
 *
 * @version $Id: ParsimonyStateStatistic.java,v 1.11 2005/07/11 14:06:25 rambaut Exp $
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 *
 */
public class ParsimonyStateStatistic extends Statistic.Abstract implements TreeStatistic {

	public static final String PARSIMONY_STATE_STATISTIC = "parsimonyStateStatistic";
	public static final String STATE = "state";
	public static final String MRCA = "mrca";

	public ParsimonyStateStatistic(String name, Tree tree, TaxonList stateTaxa, TaxonList mrcaTaxa) throws Tree.MissingTaxonException {

		super(name);
		this.tree = tree;
		this.stateLeafSet = Tree.Utils.getLeavesForTaxa(tree, stateTaxa);
		if (mrcaTaxa != null) {
			this.mrcaLeafSet = Tree.Utils.getLeavesForTaxa(tree, mrcaTaxa);
		}
	}

	public void setTree(Tree tree) { this.tree = tree; }
	public Tree getTree() { return tree; }
	
	public int getDimension() { return 1; }
	
	/** @return the parsimony tree length of the character. */
	public double getStatisticValue(int dim) {
	
		NodeRef node;
		
		if (mrcaLeafSet != null) {
			node = Tree.Utils.getCommonAncestorNode(tree, mrcaLeafSet);
		} else {
			node = tree.getRoot();
		}
		return Tree.Utils.getParsimonyState(tree, node, stateLeafSet);
	}
	
	public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {
	
		public String getParserName() { return PARSIMONY_STATE_STATISTIC; }
	
		public Object parseXMLObject(XMLObject xo) throws XMLParseException {
			
			String name;
			if (xo.hasAttribute(NAME)) {
				name = xo.getStringAttribute(NAME);
			} else {
				name = xo.getId();
			}

			Tree tree = (Tree)xo.getChild(Tree.class);
			XMLObject cxo = (XMLObject)xo.getChild(STATE);
			TaxonList stateTaxa = (TaxonList)cxo.getChild(TaxonList.class);
			
			cxo = (XMLObject)xo.getChild(MRCA);
			TaxonList mrcaTaxa = (TaxonList)cxo.getChild(TaxonList.class);	
			
			try {
				return new ParsimonyStateStatistic(name, tree, stateTaxa, mrcaTaxa);
			} catch (Tree.MissingTaxonException mte) {
				throw new XMLParseException("Taxon, " + mte + ", in " + getParserName() + "was not found in the tree.");
			}
		}
		
		//************************************************************************
		// AbstractXMLObjectParser implementation
		//************************************************************************
		
		public String getParserDescription() {
			return "A statistic that has as its value the parsimony state reconstruction of a " +
					"binary state defined by a set of taxa at a given MRCA of a tree";
		}

		public Class getReturnType() { return ParsimonyStateStatistic.class; }

		public XMLSyntaxRule[] getSyntaxRules() { return rules; }
		
		private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
			new StringAttributeRule(NAME, "A name for this statistic for the purposes of logging", true),
			new ElementRule(TreeModel.class),
			new ElementRule(STATE, 
				new XMLSyntaxRule[] { new ElementRule(Taxa.class) }),
			new ElementRule(MRCA, 
				new XMLSyntaxRule[] { new ElementRule(Taxa.class) })
		};
	};
	
	private Tree tree = null;
	private Set stateLeafSet = null;
	private Set<String> mrcaLeafSet = null;
}

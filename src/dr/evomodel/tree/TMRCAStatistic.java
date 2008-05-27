/*
 * TMRCAStatistic.java
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
 * A statistic that tracks the time of MRCA of a set of taxa
 *
 * @version $Id: TMRCAStatistic.java,v 1.21 2005/07/11 14:06:25 rambaut Exp $
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 *
 */
public class TMRCAStatistic extends Statistic.Abstract implements TreeStatistic {
	
	public static final String TMRCA_STATISTIC = "tmrcaStatistic";
	public static final String MRCA = "mrca";

	public TMRCAStatistic(String name, Tree tree, TaxonList taxa, boolean isRate) throws Tree.MissingTaxonException {
		super(name);
		this.tree = tree;
		this.leafSet = Tree.Utils.getLeavesForTaxa(tree, taxa);
        this.isRate = isRate;
	}

	public void setTree(Tree tree) { this.tree = tree; }
	public Tree getTree() { return tree; }
	
	public int getDimension() { return 1; }
	
	/** @return the height of the MRCA node. */
	public double getStatisticValue(int dim) {
	
		NodeRef node = Tree.Utils.getCommonAncestorNode(tree, leafSet);
		if (node == null) throw new RuntimeException("No node found that is MRCA of " + leafSet);
		if (isRate) {
            return tree.getNodeRate(node);
        }
        return tree.getNodeHeight(node);
	}
	
	public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {
	
		public String getParserName() { return TMRCA_STATISTIC; }
	
		public Object parseXMLObject(XMLObject xo) throws XMLParseException {
			
			String name;
			if (xo.hasAttribute(NAME)) {
				name = xo.getStringAttribute(NAME);
			} else {
				name = xo.getId();
			}
			Tree tree = (Tree)xo.getChild(Tree.class);
			TaxonList taxa = (TaxonList)xo.getElementFirstChild(MRCA);
            boolean isRate = false;
            if (xo.hasAttribute("rate")) {
                isRate = xo.getBooleanAttribute("rate");
            }

			try {
				return new TMRCAStatistic(name, tree, taxa, isRate);
			} catch (Tree.MissingTaxonException mte) {
				throw new XMLParseException("Taxon, " + mte + ", in " + getParserName() + "was not found in the tree.");
			}
		}
		
		//************************************************************************
		// AbstractXMLObjectParser implementation
		//************************************************************************

		public String getParserDescription() {
			return "A statistic that has as its value the height of the most recent common ancestor of a set of taxa in a given tree";
		}

		public Class getReturnType() { return TMRCAStatistic.class; }

		public XMLSyntaxRule[] getSyntaxRules() { return rules; }
		
		private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
			new ElementRule(TreeModel.class),
			new StringAttributeRule("name", "A name for this statistic primarily for the purposes of logging", true),
            AttributeRule.newBooleanRule("rate", true),
            new ElementRule(MRCA,
				new XMLSyntaxRule[] { new ElementRule(Taxa.class) })
		};
	};

	private Tree tree = null;
	private Set<String> leafSet = null;
    private boolean isRate;

}

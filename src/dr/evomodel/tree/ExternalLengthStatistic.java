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
import dr.evolution.util.Taxon;
import dr.inference.model.Statistic;
import dr.xml.*;

import java.util.Set;
import java.util.HashSet;
import java.util.List;
import java.util.ArrayList;

/**
 * A statistic that tracks the time of MRCA of a set of taxa
 *
 * @version $Id: TMRCAStatistic.java,v 1.21 2005/07/11 14:06:25 rambaut Exp $
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 *
 */
public class ExternalLengthStatistic extends Statistic.Abstract implements TreeStatistic {

	public static final String EXTERNAL_LENGTH_STATISTIC = "externalLengthStatistic";

	public ExternalLengthStatistic(String name, Tree tree, TaxonList taxa) throws Tree.MissingTaxonException {
		super(name);
		this.tree = tree;
        int m = taxa.getTaxonCount();
        int n = tree.getExternalNodeCount();

        for (int i = 0; i < m; i++) {

            Taxon taxon = taxa.getTaxon(i);
            NodeRef node = null;
            boolean found = false;
            for (int j = 0; j < n; j++) {

                node = tree.getExternalNode(j);
                if (tree.getNodeTaxon(node).getId().equals(taxon.getId())) {

                    found = true;
                    break;
                }
            }

            if (!found) {
                throw new Tree.MissingTaxonException(taxon);
            }

            leafSet.add(node);
        }
    }

	public void setTree(Tree tree) { this.tree = tree; }
	public Tree getTree() { return tree; }

	public int getDimension() { return leafSet.size(); }

	/** @return the height of the MRCA node. */
	public double getStatisticValue(int dim) {
        NodeRef node = leafSet.get(dim);
        return tree.getBranchLength(node);
	}

	public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

		public String getParserName() { return EXTERNAL_LENGTH_STATISTIC; }

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {

			String name;
			if (xo.hasAttribute(NAME)) {
				name = xo.getStringAttribute(NAME);
			} else {
				name = xo.getId();
			}
			Tree tree = (Tree)xo.getChild(Tree.class);
			TaxonList taxa = (TaxonList)xo.getChild(Taxa.class);

			try {
				return new ExternalLengthStatistic(name, tree, taxa);
			} catch (Tree.MissingTaxonException mte) {
				throw new XMLParseException("Taxon, " + mte + ", in " + getParserName() + "was not found in the tree.");
			}
		}

		//************************************************************************
		// AbstractXMLObjectParser implementation
		//************************************************************************

		public String getParserDescription() {
			return "A statistic that has as its value(s) the length of the external branch length(s) of a set of one or more taxa in a given tree";
		}

		public Class getReturnType() { return ExternalLengthStatistic.class; }

		public XMLSyntaxRule[] getSyntaxRules() { return rules; }

		private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
			new ElementRule(TreeModel.class),
			new StringAttributeRule("name", "A name for this statistic primarily for the purposes of logging", true),
            new ElementRule(Taxa.class)
		};
	};

	private Tree tree = null;
	private List<NodeRef> leafSet = new ArrayList<NodeRef>();

}
/*
 * MonophylyStatistic.java
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

import dr.evolution.tree.Tree;
import dr.evolution.tree.NodeRef;
import dr.evolution.util.*;
import dr.inference.model.BooleanStatistic;
import dr.xml.*;

import java.util.*;

/**
 * Tests whether 2 (possibly unresolved) trees are compatible
 * *
 * @author Andrew Rambaut
 *
 */
public class CompatibilityStatistic extends BooleanStatistic implements TreeStatistic {

	public static final String COMPATIBILITY_STATISTIC = "compatibilityStatistic";
	public static final String COMPATIBLE_WITH = "compatibleWith";

	public CompatibilityStatistic(String name, Tree tree1, Tree tree2) throws Tree.MissingTaxonException {

		super(name);
		this.tree = tree1;

		intersection = new BitSet(tree1.getExternalNodeCount());
		clades = new HashSet<BitSet>();
		getClades(tree1, tree2, tree2.getRoot(), null, clades);

		for (int i = 0; i < tree1.getTaxonCount(); i++) {
			String id = tree1.getTaxonId(i);
			if (tree2.getTaxonIndex(id) == -1) {
				throw new Tree.MissingTaxonException(tree1.getTaxon(i));
			}
		}
	}

	public void setTree(Tree tree) { this.tree = tree; }
	public Tree getTree() { return tree; }

	public int getDimension() { return 1; }

	/** @return boolean result of test. */
	public boolean getBoolean(int dim) {
		return isCompatible(tree, tree.getRoot(), null);
	}

	private boolean isCompatible(Tree tree, NodeRef node, BitSet leaves) {
		if (tree.isExternal(node)) {
			leaves.set(node.getNumber());
			return true;
		} else {

			BitSet ls = new BitSet(tree.getExternalNodeCount());

			for (int i = 0; i < tree.getChildCount(node); i++) {

				NodeRef node1 = tree.getChild(node, i);

				if (!isCompatible(tree, node1, ls)) {
					// as soon as we have an incompatibility break out...
					return false;
				}
			}

			if (leaves != null) {
				if (!clades.contains(ls)) {
					// except for the root clade...
					for (BitSet clade : clades) {
						intersection.clear();
						intersection.or(clade);
						intersection.and(ls);

						int card = intersection.cardinality();
						if (card != 0 &&
								card != ls.cardinality() &&
								card != clade.cardinality()) {
							return false;
						}
					}
				}

				leaves.or(ls);
			}
		}
		return true;
	}

	private void getClades(Tree referenceTree, Tree tree, NodeRef node, BitSet leaves, Set<BitSet> clades) {

		if (tree.isExternal(node)) {
			String taxonId = tree.getNodeTaxon(node).getId();
			for (int i = 0; i < referenceTree.getExternalNodeCount(); i++) {
				NodeRef n = referenceTree.getExternalNode(i);

				if (taxonId.equals(referenceTree.getNodeTaxon(n).getId())) {
					leaves.set(n.getNumber());
				}
			}
		} else {

			BitSet ls = new BitSet(tree.getExternalNodeCount());

			for (int i = 0; i < tree.getChildCount(node); i++) {

				NodeRef node1 = tree.getChild(node, i);

				getClades(referenceTree, tree, node1, ls, clades);
			}

			if (leaves != null) {
				// except for the root clade...
				leaves.or(ls);
				clades.add(ls);
			}

		}
	}

	public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

		public String getParserName() { return COMPATIBILITY_STATISTIC; }

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {

			String name;
			if (xo.hasAttribute(NAME)) {
				name = xo.getStringAttribute(NAME);
			} else {
				name = xo.getId();
			}

			Tree tree1 = (Tree)xo.getChild(Tree.class);

			XMLObject cxo = (XMLObject)xo.getChild(COMPATIBLE_WITH);
			Tree tree2 = (Tree)cxo.getChild(Tree.class);

			try {
				return new CompatibilityStatistic(name, tree1, tree2);
			} catch (Tree.MissingTaxonException mte) {
				throw new XMLParseException("Taxon, " + mte + ", in " + getParserName() + "was in the source tree but not the constraints tree.");
			}
		}

		//************************************************************************
		// AbstractXMLObjectParser implementation
		//************************************************************************

		public String getParserDescription() {
			return "A statistic that returns true if a pair of trees are compatible";
		}

		public Class getReturnType() { return CompatibilityStatistic.class; }

		public XMLSyntaxRule[] getSyntaxRules() { return rules; }

		private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
				new StringAttributeRule(NAME, "A name for this statistic for the purpose of logging", true),
				new ElementRule(Tree.class),
				new ElementRule(COMPATIBLE_WITH, new XMLSyntaxRule[] {
						new ElementRule(Tree.class)
				}),
				new ElementRule(Tree.class)
		};

	};

	private Tree tree;
	private final Set<BitSet> clades;
	private final BitSet intersection;

}
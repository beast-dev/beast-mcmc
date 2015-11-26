/*
 * RankProportionStatistic.java
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

package dr.app.treestat.statistics;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.TaxonList;

/**
 *
 * @version $Id: RankProportionStatistic.java,v 1.2 2005/09/28 13:50:56 rambaut Exp $
 *
 * @author Alexei Drummond
 */
public class RankProportionStatistic extends AbstractTreeSummaryStatistic {

	private RankProportionStatistic() {
		this(true);
	}

	private RankProportionStatistic(boolean proportion) {
		this.rank = 1;
		this.proportion = proportion;
	}

    public void setInteger(int value) {
        this.rank = value;
    }

    public double[] getSummaryStatistic(Tree tree) {

		double externalLength = 0.0;
		double internalLength = 0.0;
		double rankLength = 0.0;

		int externalNodeCount = tree.getExternalNodeCount();
		for (int i = 0; i < externalNodeCount; i++) {
			NodeRef node = tree.getExternalNode(i);
			NodeRef parent = tree.getParent(node);
			externalLength += tree.getNodeHeight(parent) - tree.getNodeHeight(node);
		}

		int internalNodeCount = tree.getInternalNodeCount();
		for (int i = 0; i < internalNodeCount; i++) {
			NodeRef node = tree.getInternalNode(i);
			if (!tree.isRoot(node)) {
				NodeRef parent = tree.getParent(node);
				internalLength += tree.getNodeHeight(parent) - tree.getNodeHeight(node);
			}
		}

		if (rank == 1) {
			if (!proportion) {
				return new double[] { externalLength };
			}
			return new double[] { externalLength/(internalLength+externalLength) };
		}

		for (int i = 0; i < internalNodeCount; i++) {
			NodeRef node = tree.getInternalNode(i);
			if (!tree.isRoot(node)) {
				int r = getRank(tree, node);
				if (r == rank) {
					NodeRef parent = tree.getParent(node);
					rankLength += tree.getNodeHeight(parent) - tree.getNodeHeight(node);
				}
			}
		}

		if (!proportion) return new double[] { rankLength };

		return new double[] { rankLength/(internalLength+externalLength) };
	}

	private int getRank(Tree tree, NodeRef node) {
		int childCount = tree.getChildCount(node);
		if (childCount == 0) return 1;
		int size = 0;
		for (int i = 0; i < childCount; i++) {
			size += getRank(tree, tree.getChild(node, i));
		}
		return size;
	}

	public String getSummaryStatisticName() {
		if (!proportion) return "Rank " + rank + " branch length";
		return "Rank " + rank + " branch proportion";
	}

	public String getSummaryStatisticDescription() {

		String message = "";

		if (!proportion) {
			message += "The total length of branches of rank " + rank + ". ";
		} else {
			message += "The proportion of the total length of the tree that is made up of branches of rank " + rank + ". ";
		}

		message += "A branch is rank k if it has exactly k tips below it. For example, external tips are rank 1 and internal " +
					"branches directly above cherries are rank 2.";

		return message;
	}

	public String getSummaryStatisticReference() { return FACTORY.getSummaryStatisticReference(); }
	public boolean allowsPolytomies() { return FACTORY.allowsPolytomies(); }
	public boolean allowsNonultrametricTrees() { return FACTORY.allowsNonultrametricTrees(); }
	public boolean allowsUnrootedTrees() { return FACTORY.allowsUnrootedTrees(); }
	public SummaryStatisticDescription.Category getCategory() { return FACTORY.getCategory(); }

	public static final TreeSummaryStatistic.Factory FACTORY = new TreeSummaryStatistic.Factory() {

		public TreeSummaryStatistic createStatistic() {
			return new RankProportionStatistic();
		}

		public String getSummaryStatisticName() {
			return "Rank branch proportion";
		}

		public String getSummaryStatisticDescription() {

			return "The proportion of the total length of the tree that is made up of branches of a given rank";
		}

		public String getSummaryStatisticReference() {

			return "-";
		}

		public String getValueName() { return "The rank (k):"; }

		public boolean allowsPolytomies() { return true; }

		public boolean allowsNonultrametricTrees() { return true; }

		public boolean allowsUnrootedTrees() { return true; }

		public SummaryStatisticDescription.Category getCategory() { return SummaryStatisticDescription.Category.POPULATION_GENETIC; }

        public boolean allowsWholeTree() { return true; }

        public boolean allowsCharacter() { return false; }

        public boolean allowsCharacterState() { return false; }

        public boolean allowsTaxonList() { return false; }

        public boolean allowsInteger() { return true; }

        public boolean allowsDouble() { return false; }
    };

	private int rank = 2;
	private boolean proportion = true;
}

/*
 * CollessIndex.java
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
import dr.evolution.tree.TreeUtils;

/**
 *
 * @version $Id: CollessIndex.java,v 1.2 2005/09/28 13:50:56 rambaut Exp $
 *
 * @author Alexei Drummond
 */
public class CollessIndex extends AbstractTreeSummaryStatistic {

	public CollessIndex() { }

	/**
	 * Assumes strictly bifurcating tree.
	 */
	public double[] getSummaryStatistic(Tree tree) {

		double C = 0.0;
		int n = tree.getInternalNodeCount();
		for (int i =0; i < n; i++) {

			NodeRef node = tree.getInternalNode(i);

			int r = TreeUtils.getLeafCount(tree, tree.getChild(node, 0));
			int s = TreeUtils.getLeafCount(tree, tree.getChild(node, 1));

			C += Math.abs(r-s);
		}

		n = tree.getExternalNodeCount();
		C *= 2.0 / (n * (n - 3) + 2);
		return new double[] { C };
	}

	public String getSummaryStatisticName() { return FACTORY.getSummaryStatisticName(); }
	public String getSummaryStatisticDescription() { return FACTORY.getSummaryStatisticDescription(); }
	public String getSummaryStatisticReference() { return FACTORY.getSummaryStatisticReference(); }
	public boolean allowsPolytomies() { return FACTORY.allowsPolytomies(); }
	public boolean allowsNonultrametricTrees() { return FACTORY.allowsNonultrametricTrees(); }
	public boolean allowsUnrootedTrees() { return FACTORY.allowsUnrootedTrees(); }
	public SummaryStatisticDescription.Category getCategory() { return FACTORY.getCategory(); }

	public static final TreeSummaryStatistic.Factory FACTORY = new TreeSummaryStatistic.Factory() {

		public TreeSummaryStatistic createStatistic() {
			return new CollessIndex();
		}

		public String getSummaryStatisticName() {
			return "Colless tree-imbalance";
		}

		public String getSummaryStatisticDescription() {
			return "The normalized sum of differences of number of children " +
			"in left and right subtrees over all internal nodes";
		}

		public String getSummaryStatisticReference() {
			return "Colless (1982)";
		}

		public boolean allowsPolytomies() { return false; }

		public boolean allowsNonultrametricTrees() { return true; }

		public boolean allowsUnrootedTrees() { return false; }

		public SummaryStatisticDescription.Category getCategory() { return SummaryStatisticDescription.Category.TREE_SHAPE; }
	};
}
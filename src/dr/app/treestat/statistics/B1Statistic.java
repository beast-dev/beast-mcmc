/*
 * B1Statistic.java
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

/**
 *
 * @version $Id: B1Statistic.java,v 1.2 2005/09/28 13:50:56 rambaut Exp $
 *
 * @author Alexei Drummond
 */
public class B1Statistic extends AbstractTreeSummaryStatistic {

	private B1Statistic() { }

	public double[] getSummaryStatistic(Tree tree) {

		double B1 = 0.0;
		int n = tree.getInternalNodeCount();
		for (int i =0; i < n; i++) {

			NodeRef node = tree.getInternalNode(i);
			if (!tree.isRoot(node)) {
				B1 += 1.0/getMi(tree, node);
			}
		}
		return new double[] { B1 };

	}

	/**
	 * Assumes strictly bifurcating tree
	 */
	private static int getMi(Tree tree, NodeRef node) {
		int childCount = tree.getChildCount(node);
		if (childCount == 0) return 0;
		int Mi = 0;
		for (int i =0; i < childCount; i++) {
			int mi = getMi(tree, tree.getChild(node, i));
			if (mi > Mi) Mi = mi;
		}
		Mi += 1;
		return Mi;
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
			return new B1Statistic();
		}

		public String getSummaryStatisticName() {
			return "B1";
		}

		public String getSummaryStatisticDescription() {
			return "The sum of the reciprocals of the maximum number of nodes between " +
		 	"each interior node and a tip (Mi) for all internal nodes except the " +
		 	"root.";
		}

		public String getSummaryStatisticReference() {
			return "see Kirkpatrick & Slatkin (1992)";
		}

		public boolean allowsPolytomies() { return false; }

		public boolean allowsNonultrametricTrees() { return true; }

		public boolean allowsUnrootedTrees() { return false; }

		public SummaryStatisticDescription.Category getCategory() { return SummaryStatisticDescription.Category.TREE_SHAPE; }
	};
}

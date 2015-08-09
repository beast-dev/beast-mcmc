/*
 * TreeLength.java
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
 * @version $Id: TreeLength.java,v 1.3 2005/09/28 13:50:56 rambaut Exp $
 *
 * @author Alexei Drummond
 */
public class TreeLength extends AbstractTreeSummaryStatistic {

	private TreeLength() { }

	public double[] getSummaryStatistic(Tree tree) {

		double externalLength = 0.0;
		double internalLength = 0.0;

		int externalNodeCount = tree.getExternalNodeCount();
		for (int i = 0; i < externalNodeCount; i++) {
			NodeRef node = tree.getExternalNode(i);
			externalLength += tree.getBranchLength(node);
		}

		int internalNodeCount = tree.getInternalNodeCount();
		for (int i = 0; i < internalNodeCount; i++) {
			NodeRef node = tree.getInternalNode(i);
			if (!tree.isRoot(node)) {
				internalLength += tree.getBranchLength(node);
			}
		}
		return new double[] { internalLength + externalLength };
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
			return new TreeLength();
		}

		public String getSummaryStatisticName() {
			return "Tree Length";
		}

		public String getSummaryStatisticDescription() {

			return "The sum of the branch lengths.";
		}

		public String getSummaryStatisticReference() {

			return "-";
		}

		public boolean allowsPolytomies() { return true; }

		public boolean allowsNonultrametricTrees() { return true; }

		public boolean allowsUnrootedTrees() { return true; }

		public SummaryStatisticDescription.Category getCategory() { return SummaryStatisticDescription.Category.GENERAL; }
	};
}

/*
 * TreeHeight.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.app.treestat.statistics;

import dr.evolution.tree.Tree;

/**
 *
 *
 * @author Alexei Drummond
 */
public class TreeHeight extends AbstractTreeSummaryStatistic {

	private TreeHeight() { }

	public double[] getSummaryStatistic(Tree tree) {

		return new double[] { tree.getNodeHeight(tree.getRoot()) };
	}

	public String getSummaryStatisticName() { return FACTORY.getSummaryStatisticName(); }
	public String getSummaryStatisticDescription() { return FACTORY.getSummaryStatisticDescription(); }
	public String getSummaryStatisticReference() { return FACTORY.getSummaryStatisticReference(); }
	public boolean allowsPolytomies() { return FACTORY.allowsPolytomies(); }
	public boolean allowsNonultrametricTrees() { return FACTORY.allowsNonultrametricTrees(); }
	public boolean allowsUnrootedTrees() { return FACTORY.allowsUnrootedTrees(); }
	public Category getCategory() { return FACTORY.getCategory(); }

	public static final Factory FACTORY = new Factory() {

		public TreeSummaryStatistic createStatistic() {
			return new TreeHeight();
		}

		public String getSummaryStatisticName() {
			return "Tree Height";
		}

		public String getSummaryStatisticDescription() {

			return "The height of the root of the tree.";
		}

		public String getSummaryStatisticReference() {

			return "-";
		}

		public boolean allowsPolytomies() { return true; }

		public boolean allowsNonultrametricTrees() { return true; }

		public boolean allowsUnrootedTrees() { return false; }

		public Category getCategory() { return Category.GENERAL; }
	};
}

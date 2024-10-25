/*
 * ExternalInternalRatio.java
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

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;

/**
 * @author Alexei Drummond
 *
 */
public class ExternalInternalRatio extends AbstractTreeSummaryStatistic {

	private ExternalInternalRatio() { }

	public double[] getSummaryStatistic(Tree tree) {

		double externalLength = 0.0;
		double internalLength = 0.0;

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
		return new double[] { externalLength/internalLength };
	}

	public String getSummaryStatisticName() {
		return FACTORY.getSummaryStatisticName();
	}

	public String getSummaryStatisticDescription() {
		return FACTORY.getSummaryStatisticDescription();
	}

	public String getSummaryStatisticReference() {
		return FACTORY.getSummaryStatisticReference();
	}

	public boolean allowsPolytomies() { return FACTORY.allowsPolytomies(); }

	public boolean allowsNonultrametricTrees() { return FACTORY.allowsNonultrametricTrees(); }

	public boolean allowsUnrootedTrees() { return FACTORY.allowsUnrootedTrees(); }

	public SummaryStatisticDescription.Category getCategory() { return FACTORY.getCategory(); }


	public static final TreeSummaryStatistic.Factory FACTORY = new TreeSummaryStatistic.Factory() {

		public TreeSummaryStatistic createStatistic() {
			return new ExternalInternalRatio();
		}

		public String getSummaryStatisticName() {
			return "External/Internal ratio";
		}

		public String getSummaryStatisticDescription() {

			return "The ratio of the total length of external branches to the " +
			"total length of internal branches.";
		}

		public String getSummaryStatisticReference() {
			return "-";
		}

		public boolean allowsPolytomies() { return true; }

		public boolean allowsNonultrametricTrees() { return true; }

		public boolean allowsUnrootedTrees() { return false; }

		public SummaryStatisticDescription.Category getCategory() { return SummaryStatisticDescription.Category.POPULATION_GENETIC; }

	};
}

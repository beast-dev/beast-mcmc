/*
 * IntervalKStatistic.java
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

import dr.evolution.coalescent.TreeIntervals;
import dr.evolution.tree.Tree;

/**
 * Returns the total time in the genealogy in which exactly k lineages are present.
 *
 *
 * @author Alexei Drummond
 */
public class IntervalKStatistic extends AbstractTreeSummaryStatistic {

	public IntervalKStatistic() {
		this.k = 2;
	}

    public void setInteger(int value) {
        this.k = value;
    }

    public double[] getSummaryStatistic(Tree tree) {

		TreeIntervals intervals = new TreeIntervals(tree);

		double totalTime = 0.0;
		int intervalCount = intervals.getIntervalCount();
		for (int i = 0; i < intervalCount; i++) {
			if (intervals.getLineageCount(i) == k) {
				totalTime += intervals.getInterval(i);
			}
		}
		return new double[] { totalTime };
	}

	public String getSummaryStatisticName() {
        return "TotalTime(" + k + ")";
    }

	public String getSummaryStatisticDescription() {
        return getSummaryStatisticName() + " is the total amount of time in the genealogy that " +
            "exactly " + k + " lineages exist.";
    }

	public String getSummaryStatisticReference() { return FACTORY.getSummaryStatisticReference(); }
	public boolean allowsPolytomies() { return FACTORY.allowsPolytomies(); }
	public boolean allowsNonultrametricTrees() { return FACTORY.allowsNonultrametricTrees(); }
	public boolean allowsUnrootedTrees() { return FACTORY.allowsUnrootedTrees(); }
	public SummaryStatisticDescription.Category getCategory() { return FACTORY.getCategory(); }

	public static final TreeSummaryStatistic.Factory FACTORY = new TreeSummaryStatistic.Factory() {

		public TreeSummaryStatistic createStatistic() {
			return new IntervalKStatistic();
		}

		public String getSummaryStatisticName() {
			return "TotalTime(k)";
		}

		public String getSummaryStatisticDescription() {
			return getSummaryStatisticName() + " is the total amount of time in the genealogy that " +
				"exactly k lineages exist.";
		}

		public String getSummaryStatisticReference() {
			return "-";
		}

		public String getValueName() { return "The number of lineaeges (k):"; }

		public boolean allowsPolytomies() { return true; }

		public boolean allowsNonultrametricTrees() { return true; }

		public boolean allowsUnrootedTrees() { return false; }

		public TreeSummaryStatistic.Category getCategory() { return TreeSummaryStatistic.Category.POPULATION_GENETIC; }

        public boolean allowsWholeTree() { return true; }
        public boolean allowsCharacter() { return false; }
        public boolean allowsCharacterState() { return false; }
        public boolean allowsTaxonList() { return false; }
        public boolean allowsInteger() { return true; }
        public boolean allowsDouble() { return false; }
	};

	int k = 1;
}

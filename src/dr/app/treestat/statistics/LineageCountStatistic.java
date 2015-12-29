/*
 * LineageCountStatistic.java
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

import dr.evolution.coalescent.TreeIntervals;
import dr.evolution.tree.Tree;

/**
 * Returns the total time in the genealogy in which exactly k lineages are present.
 *
 * @version $Id: IntervalKStatistic.java,v 1.2 2005/09/28 13:50:56 rambaut Exp $
 *
 * @author Alexei Drummond
 */
public class LineageCountStatistic extends AbstractTreeSummaryStatistic {

	public LineageCountStatistic() {
		this.t = 1.0;
	}

    public void setDouble(double value) {
        this.t = value;
    }

	public double[] getSummaryStatistic(Tree tree) {

		TreeIntervals intervals = new TreeIntervals(tree);

		double totalTime = 0.0;
		for (int i = 0; i < intervals.getIntervalCount(); i++) {
			totalTime += intervals.getInterval(i);
			if (totalTime > t) {
				return new double[] { intervals.getLineageCount(i) };
			}
		}
		return new double[] { 1.0 };
	}

	public String getSummaryStatisticName() {
        return "LineageCount(" + t + ")";
    }

	public String getSummaryStatisticDescription() {
        return getSummaryStatisticName() + " is the number of lineages that exists in the genealogy at " +
            "time " + t + ".";
    }

	public String getSummaryStatisticReference() { return FACTORY.getSummaryStatisticReference(); }
	public boolean allowsPolytomies() { return FACTORY.allowsPolytomies(); }
	public boolean allowsNonultrametricTrees() { return FACTORY.allowsNonultrametricTrees(); }
	public boolean allowsUnrootedTrees() { return FACTORY.allowsUnrootedTrees(); }
	public Category getCategory() { return FACTORY.getCategory(); }

	public static final Factory FACTORY = new Factory() {

		public TreeSummaryStatistic createStatistic() {
			return new LineageCountStatistic();
		}

		public String getSummaryStatisticName() {
			return "LineageCount(t)";
		}

		public String getSummaryStatisticDescription() {
			return getSummaryStatisticName() + " is the number of lineages that exists in the genealogy at " +
				"time t.";
		}

		public String getSummaryStatisticReference() {
			return "-";
		}

		public String getValueName() { return "The time (t):"; }
		public boolean allowsPolytomies() { return true; }

		public boolean allowsNonultrametricTrees() { return true; }

		public boolean allowsUnrootedTrees() { return false; }

		public Category getCategory() { return Category.POPULATION_GENETIC; }

        public boolean allowsWholeTree() { return true; }

        public boolean allowsCharacter() { return false; }

        public boolean allowsCharacterState() { return false; }

        public boolean allowsTaxonList() { return false; }

        public boolean allowsInteger() { return false; }

        public boolean allowsDouble() { return true; }
    };

	double t = 1.0;
}

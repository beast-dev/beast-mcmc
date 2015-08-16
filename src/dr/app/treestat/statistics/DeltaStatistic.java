/*
 * DeltaStatistic.java
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

import dr.evolution.tree.Tree;
import dr.math.Binomial;

import java.util.ArrayList;
import java.util.Collections;

/**
 * @author Alexei Drummond
 *
 * @version $Id: DeltaStatistic.java,v 1.1 2006/05/09 10:24:27 rambaut Exp $
 */
public class DeltaStatistic extends AbstractTreeSummaryStatistic {

    private DeltaStatistic() { }

    public double[] getSummaryStatistic(Tree tree) {

        int n = tree.getExternalNodeCount();
        double[] g = getIntervals(tree);

        double T = 0; // total branch length
        for (int j = 2; j <= n; j++) {
            T += Binomial.choose2(j) * g[j-2];
        }

        double sum = 0.0;
        for (int i = n; i > 2; i--) {
            for (int k = n; k >= i; k--) {
                sum += Binomial.choose2(k) * g[k-2];
            }
        }

        double delta = ((T / 2.0) - (sum * (1.0 / (n - 2.0)))) / (T * Math.sqrt(1.0/(12.0*(n - 2.0))));
        return new double[] { delta };
    }

    /**
     * @return the intervals in an ultrametric tree in order from root to tips.
     */
    private static double[] getIntervals(Tree tree) {

        ArrayList<Double> heights = new ArrayList<Double>();

        for (int i = 0; i < tree.getInternalNodeCount(); i++) {
            heights.add(tree.getNodeHeight(tree.getInternalNode(i)));
        }
        Collections.sort(heights, Collections.reverseOrder());

        double[] intervals = new double[heights.size()];
        for (int i = 0; i < intervals.length - 1; i++) {
            double height1 = heights.get(i);
            double height2 = heights.get(i + 1);

            intervals[i] = height1 - height2;
        }
        intervals[intervals.length - 1] = heights.get(intervals.length - 1);

        return intervals;

    }

    public String getSummaryStatisticName() { return DeltaStatistic.FACTORY.getSummaryStatisticName(); }
    public String getSummaryStatisticDescription() { return DeltaStatistic.FACTORY.getSummaryStatisticDescription(); }
    public String getSummaryStatisticReference() { return DeltaStatistic.FACTORY.getSummaryStatisticReference(); }
    public boolean allowsPolytomies() { return DeltaStatistic.FACTORY.allowsPolytomies(); }
    public boolean allowsNonultrametricTrees() { return DeltaStatistic.FACTORY.allowsNonultrametricTrees(); }
    public boolean allowsUnrootedTrees() { return DeltaStatistic.FACTORY.allowsUnrootedTrees(); }
    public Category getCategory() { return DeltaStatistic.FACTORY.getCategory(); }

    public static final Factory FACTORY = new Factory() {

        public TreeSummaryStatistic createStatistic() {
            return new DeltaStatistic();
        }

        public String getSummaryStatisticName() {
            return "Delta";
        }

        public String getSummaryStatisticDescription() {
            return "The delta-statistic is a summary of the information contained in the inter-node " +
                    "intervals of a genealogy; under the assumption of a constant-size population, " +
                    "it follows a normal distribution with mean of zero and a standard-deviation " +
                    "of one " +
                    "(Pybus and Harvey 2000). Thus, the null hypothesis that the population has a constant " +
                    "population size may be tested with 1 -  2*pnorm(abs(delta.stat(phy))) for a two-tailed test, " +
                    "or 1 - pnorm(abs(delta.stat(phy))) for a one-tailed test, both returning the corresponding P-value.";
        }

        public String getSummaryStatisticReference() {
            return "Pybus & Harvey (2000)";
        }

        public boolean allowsPolytomies() { return true; }

        public boolean allowsNonultrametricTrees() { return false; }

        public boolean allowsUnrootedTrees() { return false; }

        public Category getCategory() { return Category.POPULATION_GENETIC; }
    };
}

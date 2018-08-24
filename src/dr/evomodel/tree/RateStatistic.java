/*
 * RateStatistic.java
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

package dr.evomodel.tree;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodelxml.tree.RateStatisticParser;
import dr.inference.model.Statistic;
import dr.stats.DiscreteStatistics;

/**
 * A statistic that tracks the mean, variance and coefficent of variation of the rates.
 *
 * @author Alexei Drummond
 * @version $Id: RateStatistic.java,v 1.9 2005/07/11 14:06:25 rambaut Exp $
 */
public class RateStatistic extends TreeStatistic {

    public RateStatistic(String name, Tree tree, BranchRateModel branchRateModel, boolean external, boolean internal, String mode) {
        super(name);
        this.tree = tree;
        this.branchRateModel = branchRateModel;
        this.internal = internal;
        this.external = external;
        this.mode = mode;
    }

    public void setTree(Tree tree) {
        this.tree = tree;
    }

    public Tree getTree() {
        return tree;
    }

    public int getDimension() {
        return 1;
    }

    /**
     * @return the height of the MRCA node.
     */
    public double getStatisticValue(int dim) {

        int length = 0;
        int offset = 0;
        if (external) {
            length += tree.getExternalNodeCount();
            offset = length;
        }
        if (internal) {
            length += tree.getInternalNodeCount() - 1;
        }

        final double[] rates = new double[length];
        // need those only for mean
        final double[] branchLengths = new double[length];

        for (int i = 0; i < offset; i++) {
            NodeRef child = tree.getExternalNode(i);
            NodeRef parent = tree.getParent(child);
            branchLengths[i] = tree.getNodeHeight(parent) - tree.getNodeHeight(child);
            rates[i] = branchRateModel.getBranchRate(tree, child);
        }
        if (internal) {
            final int n = tree.getInternalNodeCount();
            int k = offset;
            for (int i = 0; i < n; i++) {
                NodeRef child = tree.getInternalNode(i);
                if (!tree.isRoot(child)) {
                    NodeRef parent = tree.getParent(child);
                    branchLengths[k] = tree.getNodeHeight(parent) - tree.getNodeHeight(child);
                    rates[k] = branchRateModel.getBranchRate(tree, child);
                    k++;
                }
            }
        }

        if (mode.equals(RateStatisticParser.MEAN)) {
            double totalWeightedRate = 0.0;
            double totalTreeLength = 0.0;
            for (int i = 0; i < rates.length; i++) {
                totalWeightedRate += rates[i] * branchLengths[i];
                totalTreeLength += branchLengths[i];
            }
            return totalWeightedRate / totalTreeLength;
        } else if (mode.equals(RateStatisticParser.VARIANCE)) {
            return DiscreteStatistics.variance(rates);
        } else if (mode.equals(RateStatisticParser.COEFFICIENT_OF_VARIATION)) {
            // don't compute mean twice
            final double mean = DiscreteStatistics.mean(rates);
            return Math.sqrt(DiscreteStatistics.variance(rates, mean)) / mean;
        }

        throw new IllegalArgumentException();
    }

    private Tree tree = null;
    private BranchRateModel branchRateModel = null;
    private boolean internal = true;
    private boolean external = true;
    private String mode = RateStatisticParser.MEAN;
}

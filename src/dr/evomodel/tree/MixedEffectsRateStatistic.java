/*
 * MixedEffectsRateStatistic.java
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

package dr.evomodel.tree;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.ArbitraryBranchRates;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodelxml.tree.MixedEffectsRateStatisticParser;
import dr.evomodelxml.tree.RateStatisticParser;
import dr.inference.model.Statistic;
import dr.stats.DiscreteStatistics;

import java.util.ArrayList;
import java.util.List;

/**
 * A statistic that tracks the mean, variance and coefficent of variation of the rates of a mixed-effects clock model.
 *
 * @author Andy Magee
 */
public class MixedEffectsRateStatistic extends TreeStatistic {

    public MixedEffectsRateStatistic(String name, Tree tree, ArbitraryBranchRates branchRateModel, boolean external, boolean internal, boolean logScale, String mode) {
        super(name);
        this.tree = tree;
        this.branchRateModel = branchRateModel;
        this.internal = internal;
        this.external = external;
        this.mode = mode;
        this.logScale = logScale;

        if (!(branchRateModel.getTransform() instanceof ArbitraryBranchRates.BranchRateTransform.LocationScaleLogNormal)) {
            throw new RuntimeException("MixedEffectsRateStatistic currently only supports LocationScaleLogNormal models.");
        }

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

    // Largely the same code as RateStatistic.getStatistic()
    private void prepareForComputation() {
        int length = 0;
        int offset = 0;
        if (external) {
            length += tree.getExternalNodeCount();
            offset = length;
        }
        if (internal) {
            length += tree.getInternalNodeCount() - 1;
        }

        rates = new double[length];
        locations = new double[length];

        ArbitraryBranchRates.BranchRateTransform.LocationScaleLogNormal branchRateTransform = ((ArbitraryBranchRates.BranchRateTransform.LocationScaleLogNormal) branchRateModel.getTransform());

        for (int i = 0; i < offset; i++) {
            NodeRef child = tree.getExternalNode(i);
            rates[i] = branchRateModel.getBranchRate(tree, child);
            locations[i] = branchRateTransform.getLocation(tree, child);
        }
        if (internal) {
            final int n = tree.getInternalNodeCount();
            int k = offset;
            for (int i = 0; i < n; i++) {
                NodeRef child = tree.getInternalNode(i);
                if (!tree.isRoot(child)) {
                    rates[k] = branchRateModel.getBranchRate(tree, child);
                    locations[k] = branchRateTransform.getLocation(tree, child);
                    k++;
                }
            }
        }
    }

    private void takeLogs() {
        for (int i = 0; i < locations.length; i++) {
            locations[i] = Math.log(locations[i]);
            rates[i] = Math.log(rates[i]);
        }
    }

    private double[] getResiduals() {
        double[] residuals = new double[locations.length];
        for (int i = 0; i < residuals.length; i++) {
            residuals[i] = rates[i] - locations[i];
        }
        return residuals;
    }

    /**
     * @return the statistic
     */
    public double getStatisticValue(int dim) {
        prepareForComputation();

        if (logScale) {
            takeLogs();
        }

        double[] residuals = getResiduals();

        if (mode.equals(MixedEffectsRateStatisticParser.MEAN)) {
            return DiscreteStatistics.mean(residuals);
        } else if (mode.equals(MixedEffectsRateStatisticParser.VARIANCE)) {
            return DiscreteStatistics.variance(residuals);
        } else if (mode.equals(MixedEffectsRateStatisticParser.PROP_EXPLAINED)) {
            return 1.0 - DiscreteStatistics.variance(residuals)/DiscreteStatistics.variance(rates);
        }

        throw new IllegalArgumentException();
    }

    private Tree tree;
    private ArbitraryBranchRates branchRateModel;
    private boolean internal;
    private boolean external;
    private String mode;
    private boolean logScale;
    private double[] rates;
    private double[] locations;
}

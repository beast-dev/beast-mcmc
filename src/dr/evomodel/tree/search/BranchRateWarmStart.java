/*
 * BranchRateWarmStart.java
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

package dr.evomodel.tree.search;

import dr.evolution.tree.NodeRef;
import dr.evomodel.branchratemodel.ArbitraryBranchRates;
import dr.evomodel.tree.TreeModel;

/**
 * Per-branch rates after a topology move (plan 7 §5.6). A branch whose length in time changes keeps its expected
 * number of substitutions, rate times time, so the move changes the tree likelihood as little as the topology allows
 * (the constant-distance idea of Zhang and Drummond 2020); each new rate stays within [0.2, 5] times the old one.
 * With {@code keep}, rates just stay with their nodes.
 * <p>
 * Rates are set on the branch-rate model's raw parameter, through the inverse of its transform: the rate of a branch
 * is transform(raw), which for TreeLime's relaxed clocks is increasing in raw (identity, a location multiplier, or a
 * location-scale log-normal map).
 *
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public final class BranchRateWarmStart {

    public enum Mode {
        PRESERVE_SUBSTITUTIONS("preserveSubstitutions"),
        KEEP("keep");

        private final String name;

        Mode(String name) {
            this.name = name;
        }

        public String getName() {
            return name;
        }

        public static Mode parse(String name) {
            for (Mode mode : values()) {
                if (mode.name.equalsIgnoreCase(name)) {
                    return mode;
                }
            }
            throw new IllegalArgumentException("Unknown rate warm start '" + name + "': use preserveSubstitutions or keep");
        }
    }

    static final double MIN_FACTOR = 0.2;
    static final double MAX_FACTOR = 5.0;

    private final TreeModel tree;
    private final ArbitraryBranchRates rates; // null with a strict clock
    private final Mode mode;

    public BranchRateWarmStart(TreeModel tree, ArbitraryBranchRates rates, Mode mode) {
        this.tree = tree;
        this.rates = rates;
        this.mode = mode;
    }

    /**
     * @return whether moves rewrite rates: per-branch rates and preserveSubstitutions
     */
    public boolean isActive() {
        return rates != null && mode == Mode.PRESERVE_SUBSTITUTIONS;
    }

    public Mode getMode() {
        return mode;
    }

    /**
     * @return the rate of the branch above the node
     */
    public double getRate(NodeRef node) {
        return rates.getBranchRate(tree, node);
    }

    /**
     * @return the rate that keeps rate times time over a branch whose time goes from before to after, within
     * [0.2, 5] times the old rate
     */
    public double preserve(double rate, double before, double after) {
        if (!(after > 0.0) || !(before > 0.0)) {
            return rate;
        }
        final double factor = Math.min(MAX_FACTOR, Math.max(MIN_FACTOR, before / after));
        return rate * factor;
    }

    /**
     * Sets the raw value whose transformed rate is the given rate.
     */
    public void setRate(NodeRef node, double rate) {
        rates.setBranchRate(tree, node, rawFor(node, rate));
    }

    // the raw value with transform(raw) = rate: exact for a transform proportional to raw, otherwise by bisection on
    // log raw (the transforms are monotone in raw)
    double rawFor(NodeRef node, double rate) {
        final ArbitraryBranchRates.BranchRateTransform transform = rates.getTransform();
        final double raw = rates.getUntransformedBranchRate(tree, node);
        final double current = transform.transform(raw, tree, node);
        if (transform instanceof ArbitraryBranchRates.BranchRateTransform.None
                || transform instanceof ArbitraryBranchRates.BranchRateTransform.MultiplyByLocation) {
            return current > 0.0 ? raw * rate / current : rate;
        }
        if (!(raw > 0.0) || !(rate > 0.0)) {
            return raw;
        }
        final boolean increasing = transform.transform(raw * 2.0, tree, node) > current;
        double lo = Math.log(raw);
        double hi = lo;
        // bracket the root in log raw
        double step = 1.0;
        for (int k = 0; k < 60 && side(transform, node, lo, rate, increasing) > 0; k++, step *= 2) {
            lo -= step;
        }
        step = 1.0;
        for (int k = 0; k < 60 && side(transform, node, hi, rate, increasing) < 0; k++, step *= 2) {
            hi += step;
        }
        for (int k = 0; k < 100 && hi - lo > 1E-13 * Math.max(1.0, Math.abs(lo)); k++) {
            final double mid = 0.5 * (lo + hi);
            if (side(transform, node, mid, rate, increasing) < 0) {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        return Math.exp(0.5 * (lo + hi));
    }

    // the sign of the (oriented) difference between the rate at exp(logRaw) and the target
    private double side(ArbitraryBranchRates.BranchRateTransform transform, NodeRef node, double logRaw, double rate,
                        boolean increasing) {
        final double difference = transform.transform(Math.exp(logRaw), tree, node) - rate;
        return increasing ? difference : -difference;
    }
}

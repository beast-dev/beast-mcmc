/*
 * LocalOptimizer.java
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
import dr.math.UnivariateFunction;
import dr.math.UnivariateMinimum;

/**
 * The search's T1 tier (plan 7 §6.1): a few node heights and branch rates around a move, optimized one coordinate at
 * a time while everything else stays fixed, so each evaluation is a local (node-event) update.
 * <p>
 * Heights are moved in bounded logit coordinates between the oldest child and the parent (log of the distance above
 * the oldest child for the root); raw rates in log coordinates. Each coordinate takes a finite-difference Newton step
 * from three evaluations (step 1E-3), clamped to 2 units; when the curvature is not positive or the step does not
 * improve, Brent's method on [y - 2, y + 2] takes over. Bounds are recomputed before each coordinate, as the
 * neighbours' heights change them. The objective never decreases: each coordinate ends at the best point evaluated.
 *
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public final class LocalOptimizer {

    /**
     * The objective, log density, at the tree's current state.
     */
    public interface Objective {
        double value();
    }

    private static final double STEP = 1E-3;
    private static final double MAX_NEWTON = 2.0;
    private static final double BRENT_RANGE = 2.0;
    private static final int BRENT_EVALUATIONS = 12;

    private final TreeModel tree;
    private final ArbitraryBranchRates rates; // null: heights only
    private final int sweeps;

    public LocalOptimizer(TreeModel tree, ArbitraryBranchRates rates, int sweeps) {
        this.tree = tree;
        this.rates = rates;
        this.sweeps = sweeps;
    }

    /**
     * @param value the objective at the current state
     * @return the objective at the end, where the state is left (never lower than value)
     */
    public double optimize(Objective objective, double value, int[] heightNodes, int[] rateNodes) {
        double f = value;
        for (int sweep = 0; sweep < sweeps; ++sweep) {
            for (int number : heightNodes) {
                f = optimize(objective, f, new HeightCoordinate(tree.getNode(number)));
            }
            if (rates != null) {
                for (int number : rateNodes) {
                    f = optimize(objective, f, new RateCoordinate(tree.getNode(number)));
                }
            }
        }
        return f;
    }

    private abstract class Coordinate {
        abstract double get();

        abstract void set(double y);

        // the model value itself (a height or a raw rate), set back exactly
        abstract double getValue();

        abstract void setValue(double value);
    }

    private final class HeightCoordinate extends Coordinate {
        private final NodeRef node;
        private final double lower;
        private final double upper;

        HeightCoordinate(NodeRef node) {
            this.node = node;
            double lower = Double.NEGATIVE_INFINITY;
            for (int c = 0; c < tree.getChildCount(node); ++c) {
                lower = Math.max(lower, tree.getNodeHeight(tree.getChild(node, c)));
            }
            this.lower = lower;
            this.upper = tree.isRoot(node) ? Double.POSITIVE_INFINITY : tree.getNodeHeight(tree.getParent(node));
        }

        double get() {
            final double h = tree.getNodeHeight(node);
            if (upper == Double.POSITIVE_INFINITY) {
                return Math.log(h - lower);
            }
            final double x = (h - lower) / (upper - lower);
            return Math.log(x / (1.0 - x));
        }

        void set(double y) {
            final double h = upper == Double.POSITIVE_INFINITY ? lower + Math.exp(y)
                    : lower + (upper - lower) / (1.0 + Math.exp(-y));
            // keep strictly inside the interval, where every branch has positive length
            if (h > lower && h < upper) {
                tree.setNodeHeight(node, h);
            }
        }

        double getValue() {
            return tree.getNodeHeight(node);
        }

        void setValue(double value) {
            tree.setNodeHeight(node, value);
        }
    }

    private final class RateCoordinate extends Coordinate {
        private final NodeRef node;

        RateCoordinate(NodeRef node) {
            this.node = node;
        }

        double get() {
            return Math.log(rates.getUntransformedBranchRate(tree, node));
        }

        void set(double y) {
            rates.setBranchRate(tree, node, Math.exp(y));
        }

        double getValue() {
            return rates.getUntransformedBranchRate(tree, node);
        }

        void setValue(double value) {
            rates.setBranchRate(tree, node, value);
        }
    }

    // minimize -F along one coordinate; returns F at the best point evaluated, where the coordinate is left
    private double optimize(final Objective objective, final double value, final Coordinate coordinate) {
        final double y0 = coordinate.get();
        if (!Double.isFinite(y0)) {
            return value;
        }
        final double[] best = {coordinate.getValue(), -value};
        final UnivariateFunction f = new UnivariateFunction() {
            public double evaluate(double y) {
                coordinate.set(y);
                double v = -objective.value();
                if (Double.isNaN(v)) {
                    v = Double.POSITIVE_INFINITY;
                }
                if (v < best[1]) {
                    best[0] = coordinate.getValue();
                    best[1] = v;
                }
                return v;
            }

            public double getLowerBound() {
                return y0 - BRENT_RANGE;
            }

            public double getUpperBound() {
                return y0 + BRENT_RANGE;
            }
        };

        final double f0 = -value;
        final double plus = f.evaluate(y0 + STEP);
        final double minus = f.evaluate(y0 - STEP);
        final double gradient = (plus - minus) / (2 * STEP);
        final double curvature = (plus - 2 * f0 + minus) / (STEP * STEP);
        boolean improved = false;
        if (curvature > 0 && Double.isFinite(gradient)) {
            final double step = Math.max(-MAX_NEWTON, Math.min(MAX_NEWTON, -gradient / curvature));
            improved = f.evaluate(y0 + step) < f0 - 1E-12 * Math.max(1.0, Math.abs(f0));
        }
        if (!improved) {
            UnivariateMinimum brent = new UnivariateMinimum();
            brent.maxFun = BRENT_EVALUATIONS;
            brent.optimize(y0, f, 1E-3, y0 - BRENT_RANGE, y0 + BRENT_RANGE);
        }
        if (coordinate.getValue() != best[0]) {
            coordinate.setValue(best[0]);
        }
        return objective.value();
    }
}

/*
 * RootToTip.java
 *
 * Copyright (C) 2002-2009 Alexei Drummond and Andrew Rambaut
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

package dr.app.pathogen;

import dr.evolution.tree.*;
import dr.evolution.util.*;
import dr.evolution.util.Date;
import dr.stats.Regression;
import dr.stats.DiscreteStatistics;
import dr.math.*;

import java.util.*;

/*
 * @author Andrew Rambaut
 */

public class TemporalRooting {

    public enum RootingFunction {
        RESIDUAL_MEAN_SQUARED("residual mean squared"),
        //        SUM_RESIDUAL_SQUARED("sum squared residuals"),
        CORRELATION("correlation"),
        R_SQUARED("R squared");

        RootingFunction(final String name) {
            this.name = name;
        }

        @Override
        public String toString() {
            return name;
        }

        private final String name;
    }

    private boolean contemporaneous = false;
    private final TaxonList taxa;
    private final Map<String, Double> dates;
    private final Map<String, Double> precisions;
    private boolean useTargetRate = false;
    private double targetRate = 0.0;
    private double dateMin;
    private double dateMax;

    private int currentRootBranch = 0;
    private int totalRootBranches = 0;

    private boolean forcePositiveRate = false;

    public TemporalRooting(TaxonList taxa) {
        this.taxa = taxa;

        dates = new HashMap<String, Double>();
        precisions = new HashMap<String, Double>();

        dateMin = Double.MAX_VALUE;
        dateMax = -Double.MAX_VALUE;

        for (int i = 0; i < taxa.getTaxonCount(); i++) {
            Taxon taxon = taxa.getTaxon(i);
            Date date = (Date)taxon.getAttribute("date");
            double d = 0.0;
            if (date != null) {
                d = date.getAbsoluteTimeValue();
                if (date.getPrecision() > 0.0) {
                    d += date.getPrecision() / 2;
                    precisions.put(taxon.getId(), date.getPrecision());
                }
            }
            if (d > dateMax) {
                dateMax = d;
            }
            if (d < dateMin) {
                dateMin = d;
            }
            dates.put(taxon.getId(), d);
        }

        if (Math.abs(dateMax - dateMin) < 1.0E-8) {
            // probably contemporaneous tips
            contemporaneous = true;
        }
    }

    public void setForcePositiveRate(boolean forcePositiveRate) {
        this.forcePositiveRate = forcePositiveRate;
    }

    public void setTargetRate(double targetRate) {
        this.targetRate = targetRate;
    }

    public boolean isContemporaneous() {
        return contemporaneous;
    }

    public double getDateRange() {
        return dateMax - dateMin;
    }

    public Tree findRoot(Tree tree, RootingFunction rootingFunction) {

        double[] dates = getTipDates(tree);
        return findGlobalRoot(tree, dates, rootingFunction, forcePositiveRate);
    }

    public Tree findLocalRoot(Tree tree, RootingFunction rootingFunction) {

        double[] dates = getTipDates(tree);
        FlexibleTree bestTree = new FlexibleTree(tree);

        findLocalRoot(bestTree, dates, rootingFunction, forcePositiveRate);

        return bestTree;
    }

    public Regression getRootToTipRegression(Tree tree) {

        if (contemporaneous) {
            throw new IllegalArgumentException("Cannot do a root to tip regression on contemporaneous tips");
        }

        double[] dates = getTipDates(tree);
        double[] distances = getRootToTipDistances(tree);

        return new Regression(dates, distances);
    }

    public Regression getNodeDensityRegression(Tree tree) {

        if (contemporaneous) {
            throw new IllegalArgumentException("Cannot do a node density regression on contemporaneous tips");
        }

        double[] dates = getTipDates(tree);
//        double[] distances = getRootToTipDistances(tree);
        double[] density = getNodeDensity(tree);

        return new Regression(dates, density);
    }

    public Regression getAncestorRootToTipRegression(Tree tree, Regression regression) {

        if (contemporaneous) {
            throw new IllegalArgumentException("Cannot do a root to tip regression on contemporaneous tips");
        }

        double[] dates = new double[tree.getExternalNodeCount()];
        double[] distances = new double[tree.getExternalNodeCount()];
        for (int i = 0; i < tree.getExternalNodeCount(); i++) {
            NodeRef tip = tree.getExternalNode(i);
            NodeRef parent = tree.getParent(tip);
            distances[i] = getRootToTipDistance(tree, parent);

            dates[i] = regression.getXIntercept() + (distances[i] / regression.getGradient());
        }

        return new Regression(dates, distances);
    }

    public double[] getRootToTipDistances(Tree tree) {

        double[] d = new double[tree.getExternalNodeCount()];
        for (int i = 0; i < tree.getExternalNodeCount(); i++) {
            NodeRef tip = tree.getExternalNode(i);
            d[i] = getRootToTipDistance(tree, tip);
        }
        return d;
    }

    public double[] getParentRootToTipDistances(Tree tree) {

        double[] d = new double[tree.getExternalNodeCount()];
        for (int i = 0; i < tree.getExternalNodeCount(); i++) {
            NodeRef tip = tree.getExternalNode(i);
            NodeRef parent = tree.getParent(tip);
            d[i] = getRootToTipDistance(tree, parent);
        }
        return d;
    }

    public double[] getRootToTipResiduals(Tree tree, Regression regression) {

        double[] r = new double[tree.getExternalNodeCount()];
        for (int i = 0; i < tree.getExternalNodeCount(); i++) {
            NodeRef tip = tree.getExternalNode(i);
            Double date = dates.get(tree.getNodeTaxon(tip).getId());
            double d = getRootToTipDistance(tree, tip);

            r[i] = regression.getResidual(date, d);
        }
        return r;
    }

    public double[] getNodeDensity(Tree tree) {

        double[] d = new double[tree.getExternalNodeCount()];
        for (int i = 0; i < tree.getExternalNodeCount(); i++) {
            NodeRef tip = tree.getExternalNode(i);
            d[i] = getNodeDensity(tree, tip);
        }
        return d;
    }

    public double[] getTipDates(Tree tree) {
        double[] d = new double[tree.getExternalNodeCount()];
        for (int i = 0; i < tree.getExternalNodeCount(); i++) {
            NodeRef tip = tree.getExternalNode(i);
            Double date = dates.get(tree.getNodeTaxon(tip).getId());
            if (date == null) {
                throw new IllegalArgumentException("Taxon, " + tree.getNodeTaxon(tip) + ", not found in taxon list");
            }
            d[i] = date;
        }
        return d;
    }

    public double[] getTipDatePrecisions(Tree tree) {
        double[] p = new double[tree.getExternalNodeCount()];
        for (int i = 0; i < tree.getExternalNodeCount(); i++) {
            NodeRef tip = tree.getExternalNode(i);
            Double precision = precisions.get(tree.getNodeTaxon(tip).getId());
            if (precision == null) {
                precision = 0.0;
            }
            p[i] = precision;
        }
        return p;
    }

    public String[] getTipLabels(Tree tree) {
        String[] labels = new String[tree.getExternalNodeCount()];
        for (int i = 0; i < tree.getExternalNodeCount(); i++) {
            NodeRef tip = tree.getExternalNode(i);
            labels[i] = tree.getNodeTaxon(tip).getId();
        }
        return labels;
    }

    private Tree findGlobalRoot(final Tree source, final double[] dates, RootingFunction rootingFunction, boolean forcePositiveRate) {

        FlexibleTree bestTree = new FlexibleTree(source);
        double minF = findLocalRoot(bestTree, dates, rootingFunction, forcePositiveRate);
        double minDiff = Double.MAX_VALUE;

        totalRootBranches = source.getNodeCount();
        for (currentRootBranch = 0; currentRootBranch < source.getNodeCount(); currentRootBranch++) {
            FlexibleTree tmpTree = new FlexibleTree(source);
            NodeRef node = tmpTree.getNode(currentRootBranch);
            if (!tmpTree.isRoot(node)) {
                double length = tmpTree.getBranchLength(node);
                tmpTree.changeRoot(node, length * 0.5, length * 0.5);

                double f = findLocalRoot(tmpTree, dates, rootingFunction, forcePositiveRate);
                if (useTargetRate) {
                    Regression r = getRootToTipRegression(tmpTree);
                    if (Math.abs(r.getGradient() - targetRate) < minDiff) {
                        minDiff = Math.abs(r.getGradient() - targetRate);
                        bestTree = tmpTree;
                    }
                } else {
                    if (f < minF) {
                        minF = f;
                        bestTree = tmpTree;
                    }
                }
            }
        }

        return bestTree;
    }

    private double findLocalRoot(final FlexibleTree tree,
                                 final double[] dates,
                                 final RootingFunction rootingFunction,
                                 final boolean forcePositiveRate) {

        if (rootingFunction == RootingFunction.RESIDUAL_MEAN_SQUARED) {
            return findAnalyticalLocalRoot(tree, dates, rootingFunction);
        }

        NodeRef node1 = tree.getChild(tree.getRoot(), 0);
        NodeRef node2 = tree.getChild(tree.getRoot(), 1);

        final double length1 = tree.getBranchLength(node1);
        final double length2 = tree.getBranchLength(node2);

        final double sumLength = length1 + length2;

        final Set<NodeRef> tipSet1 = Tree.Utils.getExternalNodes(tree, node1);
        final Set<NodeRef> tipSet2 = Tree.Utils.getExternalNodes(tree, node2);

        final double[] y = new double[tree.getExternalNodeCount()];

        UnivariateFunction f = new UnivariateFunction() {
            //        MultivariateFunction f = new MultivariateFunction() {
            public double evaluate(final double argument) {
                double l1 = argument * sumLength;

                for (NodeRef tip : tipSet1) {
                    y[tip.getNumber()] = getRootToTipDistance(tree, tip) - length1 + l1;
                }

                double l2 = (1.0 - argument) * sumLength;

                for (NodeRef tip : tipSet2) {
                    y[tip.getNumber()] = getRootToTipDistance(tree, tip) - length2 + l2;
                }

                double score;

                if (!contemporaneous) {
                    Regression r = new Regression(dates, y);

                    switch (rootingFunction) {

                        case CORRELATION:
                            score = -r.getCorrelationCoefficient();
                            break;
                        case R_SQUARED:
                            score = -r.getRSquared();
                            break;
                        case RESIDUAL_MEAN_SQUARED:
                            score = r.getResidualMeanSquared();
                            break;
                        default:
                            throw new RuntimeException("Unknown enum value");
                    }

                    if (forcePositiveRate) {
                        score = (r.getGradient() < 0.0 ? -score : score);
                    }

                } else {
                    score = DiscreteStatistics.variance(y);
                }

                return score;
            }

            public int getNumArguments() {
                return 1;
            }

            public double getLowerBound() {
                return 0;
            }

            public double getUpperBound() {
                return 1.0;
            }
        };

//        DifferentialEvolution minimum = new DifferentialEvolution(1);
//        ConjugateDirectionSearch minimum = new ConjugateDirectionSearch();
//        double[] minx = new double[] { 0.5 };
//
//        double fminx = minimum.findMinimum(f, minx);
//        double x = minx[0];

        UnivariateMinimum minimum = new UnivariateMinimum();
        double x = minimum.findMinimum(f);

        double fminx = minimum.fminx;

        double l1 = x * sumLength;
        double l2 = (1.0 - x) * sumLength;

        tree.setBranchLength(node1, l1);
        tree.setBranchLength(node2, l2);

        return fminx;
    }

    private double findAnalyticalLocalRoot(final FlexibleTree tree,
                                           final double[] t,
                                           final RootingFunction rootingFunction) {

        if (rootingFunction != RootingFunction.RESIDUAL_MEAN_SQUARED) {
            throw new UnsupportedOperationException("Analytical local root solution only for residual mean squared");
        }

        NodeRef node1 = tree.getChild(tree.getRoot(), 0);
        NodeRef node2 = tree.getChild(tree.getRoot(), 1);

        final double length1 = tree.getBranchLength(node1);
        final double length2 = tree.getBranchLength(node2);

        final double sumLength = length1 + length2;

        final Set<NodeRef> tipSet1 = Tree.Utils.getExternalNodes(tree, node1);
        final Set<NodeRef> tipSet2 = Tree.Utils.getExternalNodes(tree, node2);

        final double[] y = getRootToTipDistances(tree);

//        ## the analytical solution. Takes a vector of divergences (vd), a vector of times (vt), the children indicators (cs) and the branch length m
//        N <- length(vd)
//        n <- sum(cs)
//        d.bar <- mean(vd)
//        t.bar <- mean(vt)
//        C <- sum(vt^2) - (sum(vt)^2)/N
//        Ai <- 2*cs - (2*n-N)/N + 2*(t.bar-vt)/(C*N)*(N*sum(vt*cs) - n*sum(vt))-1
//        Bi <- vd-d.bar + (t.bar-vt)/(C*N)*(N*sum(vt*vd)- sum(vt)*sum(vd))
//        alpha <- -sum(Ai*Bi)/(m*sum(Ai^2))
//        if(bounded){# failsafe to ensure \alpha \in (0, 1)
//            alpha <- min(max(alpha, 0), 1)
//            if(n == N){alpha <- 1}
//        }
//        return(alpha)


        int N = tipSet1.size() + tipSet2.size();
        int n = tipSet1.size();

        double d_bar = DiscreteStatistics.mean(y);
        double t_bar = DiscreteStatistics.mean(t);

        final double[] cs = new double[N];
        for (NodeRef tip : tipSet1) {
            int i = tip.getNumber();
            cs[i] = 1;
        }

        double sum_tt = 0.0;
        double sum_t = 0.0;
        double sum_d = 0.0;
        double sum_td = 0.0;
        double sum_tc = 0.0;

        for (int i = 0; i < N; i++) {
            sum_tt += t[i] * t[i];
            sum_t += t[i];
            sum_d += y[i];
            sum_td += t[i] * y[i];
            sum_tc += t[i] * cs[i];
        }

        double C = sum_tt - (sum_t * sum_t / N);

        double sumAB = 0.0;
        double sumAA = 0.0;

        for (int i = 0; i < N; i++) {
            double Ai = (2 * cs[i])
                    - ((2 * n - N) / N) +
                    (2 * (t_bar - t[i]) / (C * N) * ((N * sum_tc) - (n * sum_t)) - 1);
            double Bi = (y[i] - d_bar)
                    + (t_bar - t[i]) / (C * N) * ((N * sum_td) - (sum_t * sum_d));

            sumAB += Ai * Bi;
            sumAA += Ai * Ai;
        }


        double x = - sumAB / sumLength * sumAA;
        x = Math.min(Math.max(x, 0.0), 1.0);

        double l1 = x * sumLength;
        double l2 = (1.0 - x) * sumLength;

        tree.setBranchLength(node1, l1);
        tree.setBranchLength(node2, l2);

        Regression r = new Regression(t, getRootToTipDistances(tree));

        return r.getResidualMeanSquared();
    }

    public double getRootToTipDistance(Tree tree, NodeRef node) {
        double distance = 0;
        while (node != null) {
            distance += tree.getBranchLength(node);
            node = tree.getParent(node);
        }
        return distance;
    }

    public double getNodeDensity(Tree tree, NodeRef node) {
        double density = 0;
        while (node != null) {
            density ++;
            node = tree.getParent(node);
        }
        return density;
    }

    public Tree adjustTreeToConstraints(Tree source, Map<Set<String>, double[]> cladeHeights) {

        FlexibleTree tree = new FlexibleTree(source);
        setHeightsFromDates(tree);

        adjustTreeToConstraints(tree, tree.getRoot(), null, cladeHeights);

        return tree;
    }

    public int getCurrentRootBranch() {
        return currentRootBranch;
    }

    public int getTotalRootBranches() {
        return totalRootBranches;
    }

    private double adjustTreeToConstraints(FlexibleTree tree, NodeRef node,
                                           Set<String> leaves,
                                           Map<Set<String>, double[]> cladeHeights) {

        if (!tree.isExternal(node)) {
            Set<String> l = new HashSet<String>();
            double maxChildHeight = 0.0;

            for (int i = 0; i < tree.getChildCount(node); i++) {
                NodeRef child = tree.getChild(node, i);
                double h = adjustTreeToConstraints(tree, child, l, cladeHeights);
                if (h > maxChildHeight) {
                    maxChildHeight = h;
                }
            }

            double height = tree.getNodeHeight(node);

            double lower = maxChildHeight;
            double upper = Double.POSITIVE_INFINITY;

            if (cladeHeights != null) {
                for (Set<String> clade : cladeHeights.keySet()) {
                    if (clade.equals(l)) {
                        double[] bounds = cladeHeights.get(clade);
                        lower = Math.max(bounds[0], maxChildHeight);
                        upper = bounds[1];
                    }
                }
            }

            if (lower > upper) {
                throw new IllegalArgumentException("incompatible constraints");
            }

            if (height < lower) {
                height = lower + 1E-6;
            } else if (height > upper) {
                height = (upper + lower) / 2;
            }
            tree.setNodeHeight(node, height);

            if (leaves != null) {
                leaves.addAll(l);
            }
        } else {
            leaves.add(tree.getNodeTaxon(node).getId());
        }
        return tree.getNodeHeight(node);
    }

    private void setHeightsFromDates(FlexibleTree tree) {

        dr.evolution.util.Date mostRecent = null;
        for (int i = 0; i < taxa.getTaxonCount(); i++) {
            Date date = taxa.getTaxon(i).getDate();
            if ((date != null) && (mostRecent == null || date.after(mostRecent))) {
                mostRecent = date;
            }
        }

        if (mostRecent != null) {
            TimeScale timeScale = new TimeScale(mostRecent.getUnits(), true, mostRecent.getAbsoluteTimeValue());
            double time0 = timeScale.convertTime(mostRecent.getTimeValue(), mostRecent);

            for (int i = 0; i < tree.getExternalNodeCount(); i++) {
                NodeRef tip = tree.getExternalNode(i);

                Date date = tree.getNodeTaxon(tip).getDate();
                if (date != null) {
                    tree.setNodeHeight(tip, timeScale.convertTime(date.getTimeValue(), date) - time0);
                } else {
                    tree.setNodeHeight(tip, 0.0);
                }
            }
        }
    }


}
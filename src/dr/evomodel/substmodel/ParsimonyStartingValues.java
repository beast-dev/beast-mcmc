/*
 * ParsimonyStartingValues.java
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

package dr.evomodel.substmodel;

import dr.evolution.alignment.PatternList;
import dr.evolution.datatype.DataType;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.inference.model.Parameter;
import dr.xml.Reportable;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Data-informed starting values for a CTMC with optional BSSVS indicators: asymmetric (ComplexSubstitutionModel, one
 * rate per ordered pair) or reversible (GeneralSubstitutionModel, one rate per unordered pair).
 * <p>
 * A Sankoff (unit-cost) parsimony reconstruction on the current tree gives directed transition counts n_ij.
 * <ul>
 * <li>Indicators: on for every observed transition and its reverse, plus both-way links between a hub state and every
 * strongly connected component that the observed transitions leave unreachable. The symmetric support keeps the
 * eigenvectors of Q well conditioned.</li>
 * <li>Rates: max(n, 1) / mean(n) for rates that are on, where n = n_ij (asymmetric) or n_ij + n_ji (reversible) and
 * the mean is over observed pairs; 1 (a neutral value) for rates that are off.</li>
 * <li>Clock rate: parsimony changes per site / total tree length.</li>
 * </ul>
 *
 * @author Marc A. Suchard
 */
public class ParsimonyStartingValues implements Reportable {

    private final int stateCount;
    private final double[] counts; // directed counts, [i * stateCount + j]
    private final double changesPerSite;
    private final double treeLength;
    private final boolean[] on;    // [i * stateCount + j]
    private final int observedPairs;
    private final int linkedComponents;

    private double[] indicatorValues; // last values set, for the report
    private double[] rateValues;
    private double clockValue = Double.NaN;

    public ParsimonyStartingValues(Tree tree, PatternList patterns, DataType dataType) {

        this.stateCount = dataType.getStateCount();
        this.counts = new double[stateCount * stateCount];

        double changes = 0.0;
        double sites = 0.0;
        for (int p = 0; p < patterns.getPatternCount(); ++p) {
            double weight = patterns.getPatternWeight(p);
            double[][] cost = new double[tree.getNodeCount()][];
            postOrder(tree, tree.getRoot(), patterns, p, dataType, cost);
            double[] rootCost = cost[tree.getRoot().getNumber()];
            changes += weight * preOrder(tree, tree.getRoot(), argMin(rootCost), cost, weight);
            sites += weight;
        }
        this.changesPerSite = changes / sites;

        double length = 0.0;
        for (int i = 0; i < tree.getNodeCount(); ++i) {
            NodeRef node = tree.getNode(i);
            if (!tree.isRoot(node)) {
                length += tree.getNodeHeight(tree.getParent(node)) - tree.getNodeHeight(node);
            }
        }
        this.treeLength = length;

        this.on = new boolean[stateCount * stateCount];
        int observed = 0;
        for (int i = 0; i < stateCount; ++i) {
            for (int j = 0; j < stateCount; ++j) {
                if (counts[i * stateCount + j] > 0) {
                    on[i * stateCount + j] = true;
                    ++observed;
                }
            }
        }
        this.observedPairs = observed;
        this.linkedComponents = connectComponents();

        for (int i = 0; i < stateCount; ++i) {
            for (int j = 0; j < stateCount; ++j) {
                if (on[i * stateCount + j]) {
                    on[j * stateCount + i] = true;
                }
            }
        }
    }

    private void postOrder(Tree tree, NodeRef node, PatternList patterns, int p, DataType dataType,
                           double[][] cost) {
        double[] c = new double[stateCount];
        if (tree.isExternal(node)) {
            int taxon = patterns.getTaxonIndex(tree.getNodeTaxon(node).getId());
            if (taxon < 0) {
                throw new IllegalArgumentException("Taxon " + tree.getNodeTaxon(node).getId() +
                        " is not in the pattern list");
            }
            boolean[] allowed = dataType.getStateSet(patterns.getPatternState(taxon, p));
            for (int a = 0; a < stateCount; ++a) {
                c[a] = allowed[a] ? 0.0 : Double.POSITIVE_INFINITY;
            }
        } else {
            for (int k = 0; k < tree.getChildCount(node); ++k) {
                NodeRef child = tree.getChild(node, k);
                postOrder(tree, child, patterns, p, dataType, cost);
                double[] cc = cost[child.getNumber()];
                double best = cc[argMin(cc)] + 1.0;
                for (int a = 0; a < stateCount; ++a) {
                    c[a] += Math.min(cc[a], best);
                }
            }
        }
        cost[node.getNumber()] = c;
    }

    private int preOrder(Tree tree, NodeRef node, int state, double[][] cost, double weight) {
        int changes = 0;
        for (int k = 0; k < tree.getChildCount(node); ++k) {
            NodeRef child = tree.getChild(node, k);
            double[] cc = cost[child.getNumber()];
            int best = argMin(cc);
            int childState = (cc[state] <= cc[best] + 1.0) ? state : best;
            if (childState != state) {
                counts[state * stateCount + childState] += weight;
                ++changes;
            }
            changes += preOrder(tree, child, childState, cost, weight);
        }
        return changes;
    }

    private static int argMin(double[] x) {
        int best = 0;
        for (int i = 1; i < x.length; ++i) {
            if (x[i] < x[best]) {
                best = i;
            }
        }
        return best;
    }

    /**
     * Links a hub in the largest strongly connected component to one state of every other component, both ways.
     *
     * @return the number of components that were linked
     */
    private int connectComponents() {
        List<List<Integer>> components = stronglyConnectedComponents();
        if (components.size() == 1) {
            return 0;
        }

        List<Integer> main = components.get(0);
        for (List<Integer> component : components) {
            if (component.size() > main.size()) {
                main = component;
            }
        }

        int hub = main.get(0);
        double hubDegree = -1.0;
        for (int v : main) {
            double degree = 0.0;
            for (int w = 0; w < stateCount; ++w) {
                degree += counts[v * stateCount + w] + counts[w * stateCount + v];
            }
            if (degree > hubDegree) {
                hub = v;
                hubDegree = degree;
            }
        }

        for (List<Integer> component : components) {
            if (component != main) {
                int v = component.get(0);
                on[hub * stateCount + v] = true;
                on[v * stateCount + hub] = true;
            }
        }
        return components.size() - 1;
    }

    private List<List<Integer>> stronglyConnectedComponents() { // Tarjan
        int[] index = new int[stateCount];
        int[] low = new int[stateCount];
        boolean[] onStack = new boolean[stateCount];
        Arrays.fill(index, -1);
        List<Integer> stack = new ArrayList<>();
        List<List<Integer>> components = new ArrayList<>();
        int[] counter = new int[1];
        for (int v = 0; v < stateCount; ++v) {
            if (index[v] == -1) {
                visit(v, index, low, onStack, stack, components, counter);
            }
        }
        return components;
    }

    private void visit(int v, int[] index, int[] low, boolean[] onStack, List<Integer> stack,
                       List<List<Integer>> components, int[] counter) {
        index[v] = low[v] = counter[0]++;
        stack.add(v);
        onStack[v] = true;
        for (int w = 0; w < stateCount; ++w) {
            if (on[v * stateCount + w]) {
                if (index[w] == -1) {
                    visit(w, index, low, onStack, stack, components, counter);
                    low[v] = Math.min(low[v], low[w]);
                } else if (onStack[w]) {
                    low[v] = Math.min(low[v], index[w]);
                }
            }
        }
        if (low[v] == index[v]) {
            List<Integer> component = new ArrayList<>();
            int w;
            do {
                w = stack.remove(stack.size() - 1);
                onStack[w] = false;
                component.add(w);
            } while (w != v);
            components.add(component);
        }
    }

    /**
     * Maps rate index k to (i, j): the upper triangle row-wise (GeneralSubstitutionModel, reversible), followed for
     * ComplexSubstitutionModel by the lower triangle column-wise.
     */
    private int[][] ratePairs(boolean reversible) {
        int[][] pairs = new int[reversible ? stateCount * (stateCount - 1) / 2 : stateCount * (stateCount - 1)][];
        int k = 0;
        for (int i = 0; i < stateCount; ++i) {
            for (int j = i + 1; j < stateCount; ++j) {
                pairs[k++] = new int[]{i, j};
            }
        }
        if (!reversible) {
            for (int j = 0; j < stateCount; ++j) {
                for (int i = j + 1; i < stateCount; ++i) {
                    pairs[k++] = new int[]{i, j};
                }
            }
        }
        return pairs;
    }

    public void setIndicators(Parameter indicators, boolean reversible) {
        int[][] pairs = ratePairs(reversible);
        checkDimension(indicators, pairs.length);
        indicatorValues = new double[pairs.length];
        for (int k = 0; k < pairs.length; ++k) {
            indicatorValues[k] = on[pairs[k][0] * stateCount + pairs[k][1]] ? 1.0 : 0.0;
            indicators.setParameterValueQuietly(k, indicatorValues[k]);
        }
        indicators.fireParameterChangedEvent();
    }

    /**
     * @param useSupport if true, rates outside the parsimony support are set to 1; otherwise every rate is set from
     *                   its count
     * @param reversible one rate per unordered pair, set from n_ij + n_ji
     */
    public void setRates(Parameter rates, boolean useSupport, boolean reversible) {
        int[][] pairs = ratePairs(reversible);
        checkDimension(rates, pairs.length);
        double[] pairCounts = new double[pairs.length];
        double total = 0.0;
        int observed = 0;
        for (int k = 0; k < pairs.length; ++k) {
            int i = pairs[k][0], j = pairs[k][1];
            pairCounts[k] = counts[i * stateCount + j] + (reversible ? counts[j * stateCount + i] : 0.0);
            total += pairCounts[k];
            if (pairCounts[k] > 0) {
                ++observed;
            }
        }
        double meanCount = observed > 0 ? total / observed : 1.0;
        rateValues = new double[pairs.length];
        for (int k = 0; k < pairs.length; ++k) {
            boolean inSupport = on[pairs[k][0] * stateCount + pairs[k][1]];
            rateValues[k] = (useSupport && !inSupport) ? 1.0 : Math.max(pairCounts[k], 1.0) / meanCount;
            rates.setParameterValueQuietly(k, rateValues[k]);
        }
        rates.fireParameterChangedEvent();
    }

    public void setClockRate(Parameter clockRate) {
        clockValue = getClockRate();
        clockRate.setParameterValue(0, clockValue);
    }

    private static void checkDimension(Parameter parameter, int dimension) {
        if (parameter.getDimension() != dimension) {
            throw new IllegalArgumentException("Parameter " + parameter.getId() + " has dimension " +
                    parameter.getDimension() + "; expected " + dimension);
        }
    }

    public double getClockRate() {
        return changesPerSite / treeLength;
    }

    public int getSupportSize() {
        int size = 0;
        for (boolean b : on) {
            if (b) ++size;
        }
        return size;
    }

    public String getSummary() {
        return String.format("%.1f parsimony changes per site over %d observed directed pairs; " +
                        "%d disconnected components linked to a hub; %d of %d directed pairs in the support; " +
                        "tree length %.4g, clock rate %.6g",
                changesPerSite, observedPairs, linkedComponents, getSupportSize(), stateCount * (stateCount - 1),
                treeLength, getClockRate());
    }

    @Override
    public String getReport() {
        return "Parsimony starting values: " + getSummary() + "\n" +
                "indicators: " + join(indicatorValues) + "\n" +
                "rates: " + join(rateValues) + "\n" +
                "clock: " + (Double.isNaN(clockValue) ? "not set" : String.valueOf(clockValue)) + "\n";
    }

    private static String join(double[] values) {
        if (values == null) {
            return "not set";
        }
        StringBuilder sb = new StringBuilder();
        for (double v : values) {
            sb.append(sb.length() > 0 ? " " : "").append(v);
        }
        return sb.toString();
    }
}

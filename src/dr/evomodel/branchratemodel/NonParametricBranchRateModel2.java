/*
 * TimeVaryingBranchRateModel.java
 *
 * Copyright (c) 2002-2022 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.branchratemodel;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTraitProvider;
import dr.evolution.tree.TreeUtils;
import dr.evomodel.tree.TreeModel;
import dr.math.IntegratedTransformedSplines;

import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;
import org.apache.commons.math.FunctionEvaluationException;
import org.apache.commons.math.MaxIterationsExceededException;

import java.util.Collections;
import java.util.List;
import java.util.HashMap;
import java.util.Map;

/**
 * @author Pratyusa Datta
 * @author Marc A. Suchard
 */

public class NonParametricBranchRateModel2 extends AbstractBranchRateModel
        implements DifferentiableBranchRates, Citable {

    private final Tree tree;
    private final IntegratedTransformedSplines splines;

    private boolean nodeRatesKnown;
    private boolean storedNodeRatesKnown;

    private double[] nodeRates;
    private double[] storedNodeRates;

    private static final double CACHE_PRECISION = 1e-5;


    private Map<IntegralCacheKey, Double> integralCache;
    private Map<IntegralCacheKey, Double> storedIntegralCache;

    private boolean splineParametersChanged;

    public NonParametricBranchRateModel2(String name,
                                         Tree tree,
                                         IntegratedTransformedSplines splines) {
        super(name);

        this.tree = tree;
        this.splines = splines;

        if (tree instanceof TreeModel) {
            addModel((TreeModel) tree);
        }


        addVariable(splines.getIntercept());
        addVariable(splines.getCoefficients());

        nodeRatesKnown = false;
        nodeRates = new double[tree.getNodeCount() - 1];
        integralCache = new HashMap<>();
        storedIntegralCache = null;
        splineParametersChanged = false;
    }

    @Override
    public double getBranchRate(final Tree tree, final NodeRef node) {

        assert tree == this.tree;

        if (!nodeRatesKnown) {

            NonParametricBranchRateModel2.TreeTraversal func =
                    new NonParametricBranchRateModel2.TreeTraversal.Rate(
                            nodeRates, splines, integralCache);

            calculateNodeGeneric(func);
            nodeRatesKnown = true;
        }

        return nodeRates[getParameterIndexFromNode(node)];
    }

    @Override
    public double[] updateGradientLogDensity(double[] gradientWrtBranches, double[] value, int from, int to) {

        assert from == 0;
        assert to == splines.getCoefficientDim() - 1;
        double[] gradientWrtCoefficients = new double[splines.getCoefficientDim()];

        return gradientWrtCoefficients;
    }

    interface TreeTraversal {

        void calculate(int childIndex, double start, double end);

        class Gradient implements TreeTraversal {

            Gradient(double[] gradientCoefficients, double[] gradientNodes) {
            }

            @Override
            public void calculate(int childIndex, double start, double end) {
            }
        }

        class Rate implements TreeTraversal {

            final private double[] nodeRates;
            final private IntegratedTransformedSplines approximation;
            final private Map<IntegralCacheKey, Double> integralCache;

            Rate(double[] nodeRates, IntegratedTransformedSplines approximation,
                 Map<IntegralCacheKey, Double> integralCache) {
                this.nodeRates = nodeRates;
                this.approximation = approximation;
                this.integralCache = integralCache;
            }

            @Override
            public void calculate(int childIndex, double start, double end) {
                double branchLength = end - start;
                double integral = getOrComputeIntegral(start, end);
                nodeRates[childIndex] = integral / branchLength;
            }

            private double getOrComputeIntegral(double start, double end) {

                if (!(end > start)) {
                    return 0.0;
                }

                IntegralCacheKey key = new IntegralCacheKey(roundKey(start), roundKey(end));

                Double cached = integralCache.get(key);
                if (cached != null) {
                    return cached;
                }

                double val;
                try {
                    val = approximation.getIntegral(start, end);
                } catch (FunctionEvaluationException e) {
                    throw new RuntimeException(e);
                } catch (MaxIterationsExceededException e) {
                    throw new RuntimeException(e);
                }
                integralCache.put(key, val);
                return val;
            }
        }
    }

    private void calculateNodeGeneric(NonParametricBranchRateModel2.TreeTraversal generic) {

        NodeRef root = tree.getRoot();
        double rootHeight = tree.getNodeHeight(root);

        traverseTreeByBranchGeneric(rootHeight, tree.getChild(root, 0), generic);
        traverseTreeByBranchGeneric(rootHeight, tree.getChild(root, 1), generic);
    }

    private void traverseTreeByBranchGeneric(double currentHeight, NodeRef child,
                                             NonParametricBranchRateModel2.TreeTraversal generic) {

        final double childHeight = tree.getNodeHeight(child);
        final int childIndex = getParameterIndexFromNode(child);
        double start = childHeight;
        double end = currentHeight;

        if (end > start) {
            generic.calculate(childIndex, start, end);
        }

        if (!tree.isExternal(child)) {
            traverseTreeByBranchGeneric(childHeight, tree.getChild(child, 0), generic);
            traverseTreeByBranchGeneric(childHeight, tree.getChild(child, 1), generic);
        }
    }

    private static final boolean ROUND = true;

    private static double roundKey(double x) {
        if (ROUND) {
            if (Double.isNaN(x) || Double.isInfinite(x)) {
                return x;
            }
            return Math.round(x / CACHE_PRECISION) * CACHE_PRECISION;
        } else {
            return x;
        }
    }

    private final static class IntegralCacheKey {

        private final double lower;
        private final double upper;

        IntegralCacheKey(double lower, double upper) {
            this.lower = lower;
            this.upper = upper;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (!(o instanceof IntegralCacheKey)) return false;
            IntegralCacheKey other = (IntegralCacheKey) o;
            return Double.compare(lower, other.lower) == 0 &&
                    Double.compare(upper, other.upper) == 0;
        }

        @Override
        public int hashCode() {
            long bits1 = Double.doubleToLongBits(lower);
            long bits2 = Double.doubleToLongBits(upper);
            long combined = bits1 ^ (bits2 << 1);
            return (int)(combined ^ (combined >>> 32));
        }
    }

    @Override
    public double getBranchRateDifferential(final Tree tree, final NodeRef node) {
        return 1.0;
    }

    @Override
    public double getBranchRateSecondDifferential(Tree tree, NodeRef node) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public void handleModelChangedEvent(Model model, Object object, int index) {
        nodeRatesKnown = false;
        fireModelChanged();

        if (model != tree) {
            throw new IllegalArgumentException("How did we get here?");
        }
    }

    @Override
    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        nodeRatesKnown = false;

        if (variable == splines.getIntercept() || variable == splines.getCoefficients()) {
            integralCache = new HashMap<>();
            splineParametersChanged = true;
        }

        fireModelChanged();
    }

    @Override
    protected void storeState() {

        if (storedNodeRates == null) {
            storedNodeRates = new double[nodeRates.length];
        }

        System.arraycopy(nodeRates, 0, storedNodeRates, 0, nodeRates.length);
        storedNodeRatesKnown = nodeRatesKnown;

        storedIntegralCache = integralCache;
        splineParametersChanged = false;
    }

    @Override
    protected void restoreState() {
        double[] tmp = nodeRates;
        nodeRates = storedNodeRates;
        storedNodeRates = tmp;
        nodeRatesKnown = storedNodeRatesKnown;

        if (splineParametersChanged) {
            integralCache = storedIntegralCache;
            splineParametersChanged = false;
        }
    }

    @Override
    protected void acceptState() {
        storedIntegralCache = null;
        splineParametersChanged = false;
    }

    @Override
    public Parameter getRateParameter() {
        return splines.getCoefficients();
    }

    @Override
    public int getParameterIndexFromNode(NodeRef node) {
        int nodeNumber = node.getNumber();
        if (nodeNumber > tree.getRoot().getNumber()) {
            --nodeNumber;
        }
        return nodeNumber;
    }

    @Override
    public ArbitraryBranchRates.BranchRateTransform getTransform() {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public double[] updateDiagonalHessianLogDensity(double[] diagonalHessian, double[] gradient, double[] value,
                                                    int from, int to) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public Citation.Category getCategory() {
        return Citation.Category.MOLECULAR_CLOCK;
    }

    @Override
    public String getDescription() {
        return "Nonparametric branch rate model using splines";
    }

    @Override
    public List<Citation> getCitations() {
        return Collections.singletonList(
                new Citation(
                        new Author[]{
                                new Author("P", "Datta"),
                                new Author("P", "Lemey"),
                                new Author("MA", "Suchard"),
                        },
                        Citation.Status.IN_PREPARATION
                )
        );
    }

    public String toString() {
        TreeTraitProvider[] treeTraitProviders = {this};
        return TreeUtils.newick(tree, treeTraitProviders);
    }
}
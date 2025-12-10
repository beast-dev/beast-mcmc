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

import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTraitProvider;
import dr.evolution.tree.TreeUtils;
import dr.evomodel.bigfasttree.BigFastTreeIntervals;
import dr.evomodel.tree.TreeModel;
import dr.math.IntegratedSquaredSplines;

import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.math.matrixAlgebra.RobustEigenDecomposition;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;
import org.apache.commons.math.FunctionEvaluationException;
import org.apache.commons.math.MaxIterationsExceededException;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.HashMap;
import java.util.Map;
import java.util.Objects;

/**
 * @author Pratyusa Datta
 * @author Marc A. Suchard
 */

// This branch rate model uses quadratic splines and gives accurate rates and gradients
// Also being implemented as the Integrable FunctionalForm in timeVaryingRates
// Current Integrable version in timeVaryingRates requires debugging
// This branch rate model will no longer be needed once merged with timeVaryingRates

public class NonParametricBranchRateModel2 extends AbstractBranchRateModel
        implements DifferentiableBranchRates, Citable {



    private final Tree tree;
    private final Parameter coefficients;
    private final double[] knots;
    private final int degree;

    private boolean nodeRatesKnown;
    private boolean storedNodeRatesKnown;
    private boolean coefficientsChanged;
    private boolean restoreCache;

    private double[] nodeRates;
    private double[] storedNodeRates;
    private Map<IntegralCacheKey, Double> integralCache;
    private Map<IntegralCacheKey, Double> savedIntegralCache;
    private static final double CACHE_PRECISION = 1e-5;

    public NonParametricBranchRateModel2(String name,
                                        Tree tree,
                                        Parameter coefficients,
                                        double[] knots,
                                        int degree) {
        super(name);

        this.tree = tree;
        this.coefficients = coefficients;
        this.degree = degree;
        this.knots = knots;

        if (tree instanceof TreeModel) {
            addModel((TreeModel) tree);
        }

        addVariable(coefficients);



        nodeRatesKnown = false;
        coefficientsChanged = true;
        nodeRates = new double[tree.getNodeCount() - 1];
        integralCache = new HashMap<>();
        restoreCache = false;

    }


    @Override
    public double getBranchRate(final Tree tree, final NodeRef node) {

        assert tree == this.tree;

        if (coefficientsChanged) {
            savedIntegralCache = integralCache;
            integralCache = new HashMap<>();

            coefficientsChanged = false;
            restoreCache = true;
        }

        if (!nodeRatesKnown) {

            NonParametricBranchRateModel2.TreeTraversal func =
                    new NonParametricBranchRateModel2.TreeTraversal.Rate(
                            nodeRates, new IntegratedSquaredSplines(
                                    coefficients.getParameterValues(), knots, degree),
                            integralCache);

            calculateNodeGeneric(func);
            nodeRatesKnown = true;
        }

        return nodeRates[getParameterIndexFromNode(node)];
    }



    @Override
    public double[] updateGradientLogDensity(double[] gradientWrtBranches, double[] value, int from, int to) {

        assert from == 0;
        assert to == coefficients.getDimension() - 1;
        double[] gradientWrtCoefficients = new double[coefficients.getDimension()];

        Arrays.fill(gradientWrtCoefficients, 0);

        NonParametricBranchRateModel2.TreeTraversal func =
                new NonParametricBranchRateModel2.TreeTraversal.Gradient(gradientWrtCoefficients,
                        gradientWrtBranches, new IntegratedSquaredSplines(
                        coefficients.getParameterValues(), knots, degree));
        calculateNodeGeneric(func);


        return gradientWrtCoefficients;

    }
    interface TreeTraversal {

        void calculate(int childIndex, double start, double end, IntegratedSquaredSplines approximation,
                       double[] coefficients);

        double[] getCoefficients();

        class Gradient implements TreeTraversal {

            final private double[] gradientCoefficients;
            final private double[] gradientNodes;
            final private IntegratedSquaredSplines approximation;


            Gradient(double[] gradientCoefficients, double[] gradientNodes,
                     IntegratedSquaredSplines approximation) {

                this.gradientCoefficients = gradientCoefficients;
                this.gradientNodes = gradientNodes;
                this.approximation = approximation;
            }

            @Override
            public double[] getCoefficients() { return gradientCoefficients; }

            @Override
            public void calculate(int childIndex, double start, double end,
                                  IntegratedSquaredSplines dead_approximation, double[] coefficients) {

                int dim = getCoefficients().length;
                double branchLength = end - start;
                double[] gradientBranchWrtCoeff = approximation.getGradient(start, end);
                for (int i = 0; i < dim; i++) {
                    gradientCoefficients[i] += (gradientNodes[childIndex] * gradientBranchWrtCoeff[i])/branchLength;
                }
            }
        }

        class Rate implements TreeTraversal {

            final private double[] nodeRates;
            final private IntegratedSquaredSplines approximation;
            final private Map<IntegralCacheKey, Double> integralCache;

            Rate(double[] nodeRates, IntegratedSquaredSplines approximation, Map<IntegralCacheKey, Double> integralCache) {
                this.nodeRates = nodeRates;
                this.approximation = approximation;
                this.integralCache = integralCache;
            }

            @Override
            public double[] getCoefficients() {
                return null;
            }

            @Override
            public void calculate(int childIndex, double start, double end,
                                  IntegratedSquaredSplines dead_approximation, double[] dead_coefficients) {


                double branchLength = end - start;

                if (!USE_CACHE) {
                    double integral = approximation.getIntegral(start, end);
                    nodeRates[childIndex] = integral / branchLength;
                } else {
//                    try {
                        double integral = getOrComputeIntegral(null, start, end, approximation, integralCache);
                        nodeRates[childIndex] = integral / branchLength;
//                    } catch (FunctionEvaluationException e) {
//                        throw new RuntimeException(e);
//                    } catch (MaxIterationsExceededException e) {
//                        throw new RuntimeException(e);
//                    }
                }

            }

            private double getOrComputeIntegral(double[] dead_coeffValues, double start, double end,
                                                IntegratedSquaredSplines approximation
//    )
                    , Map<IntegralCacheKey, Double> integralCache)
//            throws FunctionEvaluationException, MaxIterationsExceededException
            {

                if (!(end > start)) {
                    return 0.0;
                }

                IntegralCacheKey key = new IntegralCacheKey(approximation.coefficient, roundKey(start), roundKey(end));

                Double cached = integralCache.get(key);
                if (cached != null) {
                    return cached;
                }

                double val = approximation.getIntegral(start, end);
                integralCache.put(key, val);
                return val;
            }

            final private static boolean USE_CACHE = true;
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
//        double[] coefficientValues = new double[coefficients.getDimension()];
        double[] coefficientValues = generic.getCoefficients();
        double start = childHeight;
        double end = currentHeight;

//        for (int i = 0; i < coefficients.getDimension(); ++i) {
//            coefficientValues[i] = coefficients.getParameterValue(i);
//        }

//        IntegratedSquaredSplines approximation = new IntegratedSquaredSplines(coefficientValues, knots, degree);

        if (end > start) {

            generic.calculate(childIndex, start, end, null, null);
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
        private final double[] coeffCopy;
        private final double lower;
        private final double upper;
        private final int coefHash;

        IntegralCacheKey(double[] dead_coefficients, double lower, double upper) {
//            this.coeffCopy = Arrays.copyOf(coefficients, coefficients.length);
            this.coeffCopy = null;
            this.lower = lower;
            this.upper = upper;
//            this.coefHash = Arrays.hashCode(this.coeffCopy);
            this.coefHash = -1;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (!(o instanceof IntegralCacheKey)) return false;
            IntegralCacheKey other = (IntegralCacheKey) o;
            return Double.compare(lower, other.lower) == 0 &&
                    Double.compare(upper, other.upper) == 0 //&&
                    //Arrays.equals(coeffCopy, other.coeffCopy)
                    ;
        }

        private static final boolean NEW_HASH = true;

        @Override
        public int hashCode() {
            if (NEW_HASH) {
                long bits1 = Double.doubleToLongBits(lower);
                long bits2 = Double.doubleToLongBits(upper);
                long combined = bits1 ^ (bits2 << 1);
                return (int)(combined ^ (combined >>> 32));
            } else {
                return Objects.hash(lower, upper);
            }
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

        if (variable == coefficients) {
            coefficientsChanged = true;
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
    }

    @Override
    protected void restoreState() {
        double[] tmp = nodeRates;
        nodeRates = storedNodeRates;
        storedNodeRates = tmp;

        nodeRatesKnown = storedNodeRatesKnown;

        if (restoreCache) {
            integralCache = savedIntegralCache;
            restoreCache = false;
        }
        restoreCache = false;
        coefficientsChanged = false;
    }

    @Override
    protected void acceptState() {
        restoreCache = false;
    }

    @Override
    public Parameter getRateParameter() {
        return coefficients;
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
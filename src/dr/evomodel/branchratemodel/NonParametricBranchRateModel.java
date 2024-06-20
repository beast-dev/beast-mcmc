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
import dr.math.IntegratedSquaredGPApproximation;


import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.math.matrixAlgebra.RobustEigenDecomposition;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * @author Pratyusa Datta
 * @author Marc A. Suchard
 */
// TO DO: Remove code duplication
public class NonParametricBranchRateModel extends AbstractBranchRateModel
        implements DifferentiableBranchRates, Citable {



    private final Tree tree;
    private final Parameter coefficients;
    private final Parameter boundary;
    private final Parameter scaleFactor;
    private final Parameter lengthScale;

    private final BigFastTreeIntervals intervals;
    private boolean nodeRatesKnown;
    private boolean storedNodeRatesKnown;

    private double[] nodeRates;
    private double[] gradientwrtCoefficients;
    private double[] storedNodeRates;
    private double[] y;

    public NonParametricBranchRateModel(String name,
                                        Tree tree,
                                        Parameter coefficients,
                                        Parameter boundary,
                                        Parameter scaleFactor,
                                        Parameter lengthScale) {
        super(name);

        this.tree = tree;
        this.coefficients = coefficients;
        this.boundary = boundary;
        this.scaleFactor = scaleFactor;
        this.lengthScale = lengthScale;

        if (tree instanceof TreeModel) {
            addModel((TreeModel) tree);
        }

        addVariable(coefficients);
        addVariable(boundary);
        addVariable(scaleFactor);
        addVariable(lengthScale);

        intervals = new BigFastTreeIntervals((TreeModel) tree);
        addModel(intervals);



        nodeRatesKnown = false;
        nodeRates = new double[tree.getNodeCount() - 1];
        gradientwrtCoefficients = new double[coefficients.getDimension()];
        y = new double[tree.getNodeCount() - 1];
    }


    private void calculateNodeRates() {

        NodeRef root = tree.getRoot();
        double rootHeight = tree.getNodeHeight(root);

        getNodeRateByBranch(rootHeight, tree.getChild(root, 0));
        getNodeRateByBranch(rootHeight, tree.getChild(root, 1));
    }

    private void getNodeRateByBranch(double currentHeight, NodeRef child) {

        final double childHeight = tree.getNodeHeight(child);
        final double branchLength = currentHeight - childHeight;
        final int childIndex = getParameterIndexFromNode(child);
        double[] coefficientValues = new double[coefficients.getDimension()];

        for (int i = 0; i < coefficients.getDimension(); ++i) {
            coefficientValues[i] = coefficients.getParameterValue(i);
        }

        IntegratedSquaredGPApproximation approximation = new IntegratedSquaredGPApproximation(coefficientValues, boundary.getParameterValue(0),
                                            scaleFactor.getParameterValue(0), lengthScale.getParameterValue(0));
        nodeRates[childIndex] = approximation.getIntegral(childHeight, currentHeight)/branchLength;

        if (!tree.isExternal(child)) {
            getNodeRateByBranch(childHeight, tree.getChild(child, 0));
            getNodeRateByBranch(childHeight, tree.getChild(child, 1));
        }
    }

    @Override
    public double getBranchRate(final Tree tree, final NodeRef node) {

        assert tree == this.tree;

        if (!nodeRatesKnown) {
            calculateNodeRates();
            nodeRatesKnown = true;
        }

        return nodeRates[getParameterIndexFromNode(node)];
    }



//    private double[] getGradientBranchWrtCoefficient(int i) {
//        NodeRef root = tree.getRoot();
//        double rootHeight = tree.getNodeHeight(root);
//
//        gradientCoefficientByBranch(rootHeight, tree.getChild(root, 0), i);
//        gradientCoefficientByBranch(rootHeight, tree.getChild(root, 1), i);
//
//        return y;
//    }



// TO DO: gradient with respect to rates
    @Override
    public double[] updateGradientLogDensity(double[] gradientWrtBranches, double[] value, int from, int to) {

        assert from == 0;
        assert to == coefficients.getDimension() - 1;

        int m = coefficients.getDimension();
        NodeRef root = tree.getRoot();
        NodeRef child1 = tree.getChild(root, 0);
        NodeRef child2 = tree.getChild(root, 1);
        double rootHeight = tree.getNodeHeight(root);
        Arrays.fill(gradientwrtCoefficients, 0);



        return gradientwrtCoefficients;

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
    }

    @Override
    protected void acceptState() { }

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
        return "Time-varying branch rate model";
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

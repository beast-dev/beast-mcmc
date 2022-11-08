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
import dr.evomodelxml.branchratemodel.TimeVaryingBranchRateModelParser;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;

import java.util.Collections;
import java.util.List;

/**
 * @author Pratyusa Datta
 * @author Marc A. Suchard
 */

public class TimeVaryingBranchRateModel extends AbstractBranchRateModel implements DifferentiableBranchRates, Citable {

    private final Tree tree;
    private final Parameter rates;
    private final Parameter gridPoints;

    private boolean nodeRatesKnown;
    private boolean storedNodeRatesKnown;

    private double[] nodeRates;
    private double[] storedNodeRates;

    private final double[] times;

    public TimeVaryingBranchRateModel(Tree tree,
                                      Parameter rates,
                                      Parameter gridPoints) {

        super(TimeVaryingBranchRateModelParser.PARSER_NAME);

        this.tree = tree;
        this.rates = rates;
        this.gridPoints = gridPoints;

        if (tree instanceof TreeModel) {
            addModel((TreeModel) tree);
        }

        addVariable(rates);
        addVariable(gridPoints);

        nodeRates = new double[tree.getNodeCount()];
        storedNodeRates = new double[tree.getNodeCount()];

        times = computeTimes();

        nodeRatesKnown = false;
    }

    @Override
    public double getBranchRate(final Tree tree, final NodeRef node) {

        assert tree == this.tree;

        if (!nodeRatesKnown) { // lazy evaluation
            calculateNodeRates();
            nodeRatesKnown = true;
        }

        return nodeRates[getParameterIndexFromNode(node)];
    }

    @Override
    public double[] updateGradientLogDensity(double[] gradient, double[] value, int from, int to) {
        // TODO
        return new double[rates.getDimension()];
    }

    @Override
    public double getBranchRateDifferential(final Tree tree, final NodeRef node) {
        return 1.0; // TODO
    }

    @Override
    public double getBranchRateSecondDifferential(Tree tree, NodeRef node) {
        throw new RuntimeException("Not yet implemented");
    }

    private void calculateNodeRates() {

        NodeRef root = tree.getRoot();
        double rootHeight = tree.getNodeHeight(root);

        int epochIndex = times.length - 1;
        while (times[epochIndex] >= rootHeight) {
            --epochIndex;
        }

        traverseTreeByBranch(rootHeight, tree.getChild(root, 0), epochIndex);
        traverseTreeByBranch(rootHeight, tree.getChild(root, 1), epochIndex);
    }

    public String toString() {
        TreeTraitProvider[] treeTraitProviders = {this};
        return TreeUtils.newick(tree, treeTraitProviders);
    }

    private void traverseTreeByBranch(double parentHeight, NodeRef child, int epochIndex) {

        // TODO needs testing / debugging

        final double childHeight = tree.getNodeHeight(child);

        double currentHeight = parentHeight;
        double branchRateNumerator = 0.0;
        double branchRateDenominator = 0.0;

        final double weightedRate;
        if (currentHeight > childHeight) {

            while (times[epochIndex] > childHeight) {
                double timeLength = currentHeight - times[epochIndex];
                double rate = rates.getParameterValue(epochIndex);

                branchRateNumerator += rate * timeLength;
                branchRateDenominator += timeLength;
                currentHeight = times[epochIndex];

                --epochIndex;
            }

            double timeLength = currentHeight - childHeight;
            double rate = rates.getParameterValue(epochIndex);

            branchRateNumerator += rate * timeLength;
            branchRateDenominator += timeLength;

            weightedRate = branchRateNumerator / branchRateDenominator;
        } else {
            weightedRate = times[epochIndex];
        }

        nodeRates[getParameterIndexFromNode(child)] = weightedRate;

        if (!tree.isExternal(child)) {
            traverseTreeByBranch(childHeight, tree.getChild(child, 0), epochIndex);
            traverseTreeByBranch(childHeight, tree.getChild(child, 1), epochIndex);
        }
    }

    private double[] computeTimes() {
        double[] times = new double[rates.getDimension()];
        System.arraycopy(gridPoints.getParameterValues(), 0, times, 1, gridPoints.getDimension());
        return times;
    }

    @Override
    public void handleModelChangedEvent(Model model, Object object, int index) {
        nodeRatesKnown = false;
        fireModelChanged();
    }

    @Override
    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {

        if (variable == gridPoints) {
            throw new RuntimeException("Not yet implemented");
        }

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
        return rates;
    }

    @Override
    public int getParameterIndexFromNode(NodeRef node) {
        return node.getNumber(); // TODO Unsure if this is correct
    }

    @Override
    public ArbitraryBranchRates.BranchRateTransform getTransform() {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public double[] updateDiagonalHessianLogDensity(double[] diagonalHessian, double[] gradient, double[] value, int from, int to) {
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
                                new Author("MA", "Suchard"),
                        },
                        Citation.Status.IN_PREPARATION
                )
        );
    }
}
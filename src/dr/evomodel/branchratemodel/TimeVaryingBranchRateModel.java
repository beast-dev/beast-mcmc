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
    private final FunctionalForm functionalForm;

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
        functionalForm = new FunctionalForm.PiecewiseConstant(rates);

        nodeRatesKnown = false;
    }

    @Override
    public double getBranchRate(final Tree tree, final NodeRef node) {

        assert tree == this.tree;

        if (!nodeRatesKnown) { // lazy evaluation
            Traversal func = new Traversal.Rates(nodeRates, functionalForm);
            calculateNodeGeneric(func);
            nodeRatesKnown = true;
        }

        return nodeRates[getParameterIndexFromNode(node)];
    }

    @Override
    public double[] updateGradientLogDensity(double[] gradientWrtBranches, double[] value, int from, int to) {

        assert from == 0;
        assert to == rates.getDimension() - 1;

        double[] gradientWrtRates = new double[rates.getDimension()];

        Traversal func = new Traversal.Gradient(gradientWrtRates, gradientWrtBranches, functionalForm);
        calculateNodeGeneric(func);

        return gradientWrtRates; 
    }

    @Override
    public double getBranchRateDifferential(final Tree tree, final NodeRef node) {
        return 1.0;
    }

    @Override
    public double getBranchRateSecondDifferential(Tree tree, NodeRef node) {
        throw new RuntimeException("Not yet implemented");
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

    interface FunctionalForm {

        void reset();

        void incrementRate(int epochIndex, double startTime, double endTime);

        double gradientWeight(int epochIndex, double startTime, double endTime, double branchLength);

        double getRateParameter(int epochIndex);

        double rateNumerator();

        abstract class Base implements FunctionalForm {

            final Parameter parameter;

            Base(Parameter parameter) {
                this.parameter = parameter;
            }
        }

        class PiecewiseConstant extends Base {

            private double branchRateNumerator;

            PiecewiseConstant(Parameter parameter) {
                super(parameter);
            }

            @Override
            public void reset() {
                branchRateNumerator = 0.0;
            }

            @Override
            public double getRateParameter(int epochIndex) {
                return parameter.getParameterValue(epochIndex);
            }

            @Override
            public void incrementRate(int epochIndex, double startTime, double endTime) {
                double timeLength = startTime - endTime;
                branchRateNumerator += getRateParameter(epochIndex) * timeLength;
            }

            @Override
            public double rateNumerator() {
                return branchRateNumerator;
            }

            @Override
            public double gradientWeight(int epochIndex, double startTIme, double endTime, double branchLength) {
                double timeLength = startTIme - endTime;
                return timeLength / branchLength;
            }
        }

        @SuppressWarnings("unused")
        class PiecewiseLogConstant extends PiecewiseConstant {

            PiecewiseLogConstant(Parameter parameter) {
                super(parameter);
            }

            @Override
            public double getRateParameter(int epochIndex) {
                return Math.exp(super.getRateParameter(epochIndex));
            }

            @Override
            public double gradientWeight(int epochIndex, double startTime, double endTime, double branchLength) {
                return super.gradientWeight(epochIndex, startTime, endTime, branchLength) *
                        getRateParameter(epochIndex);  // TODO Currently untested
            }
        }

        @SuppressWarnings("unused")
        abstract class PiecewiseLinear implements FunctionalForm { }

        @SuppressWarnings("unused")
        abstract class Integrable extends Base {

            Integrable(Parameter parameter) {
                super(parameter);
            }
        }
    }

    interface Traversal {

        void reset();

        void increment(int epochIndex, int childIndex, double startTime, double endTime, double branchLength);

        void store(int epochIndex, int nodeIndex, double branchLength);

        abstract class AbstractTraversal implements Traversal {

            final FunctionalForm functionalForm;

            AbstractTraversal(FunctionalForm functionalForm) {
                this.functionalForm = functionalForm;
            }

            @Override
            public void reset() {
                functionalForm.reset();
            }
        }

        class Gradient extends AbstractTraversal {

            final private double[] gradientEpochs;
            final private double[] gradientNodes;

            Gradient(double[] gradientEpochs, double[] gradientNodes, FunctionalForm functionalForm) {
                super(functionalForm);
                this.gradientEpochs = gradientEpochs;
                this.gradientNodes = gradientNodes;
            }

            @Override
            public void increment(int epochIndex, int childIndex,
                                  double startTime, double endTime, double branchLength) {
                gradientEpochs[epochIndex] += gradientNodes[childIndex] *
                        functionalForm.gradientWeight(epochIndex, startTime, endTime, branchLength);
            }

            @Override
            public void store(int epochIndex, int nodeIndex, double rate) { }
        }

        class Rates extends AbstractTraversal {

            final private double[] nodeRates;

            Rates(double[] nodeRates, FunctionalForm functionalForm) {
                super(functionalForm);
                this.nodeRates = nodeRates;
            }

            @Override
            public void increment(int epochIndex, int childIndex,
                                  double startTime, double endTime, double branchLength) {
                functionalForm.incrementRate(epochIndex, startTime, endTime);
            }

            @Override
            public void store(int epochIndex, int nodeIndex, double branchLength) {
                nodeRates[nodeIndex] = functionalForm.rateNumerator() / branchLength;
            }
        }
    }

    private void calculateNodeGeneric(Traversal generic) {

        NodeRef root = tree.getRoot();
        double rootHeight = tree.getNodeHeight(root);

        int epochIndex = times.length - 1;
        while (times[epochIndex] >= rootHeight) {
            --epochIndex;
        }

        traverseTreeByBranchGeneric(rootHeight, tree.getChild(root, 0), epochIndex, generic);
        traverseTreeByBranchGeneric(rootHeight, tree.getChild(root, 1), epochIndex, generic);
    }

    private void traverseTreeByBranchGeneric(double currentHeight, NodeRef child, int epochIndex,
                                             Traversal generic) {

        final double childHeight = tree.getNodeHeight(child);
        final double branchLength = currentHeight - childHeight;
        final int childIndex =  getParameterIndexFromNode(child);

        generic.reset();

        if (currentHeight > childHeight) {

            while (times[epochIndex] > childHeight) {
                generic.increment(epochIndex, childIndex, currentHeight, times[epochIndex], branchLength);
                currentHeight = times[epochIndex];

                --epochIndex;
            }

            generic.increment(epochIndex, childIndex, currentHeight, childHeight, branchLength);
        }

        generic.store(epochIndex, childIndex, branchLength);

        if (!tree.isExternal(child)) {
            traverseTreeByBranchGeneric(childHeight, tree.getChild(child, 0), epochIndex, generic);
            traverseTreeByBranchGeneric(childHeight, tree.getChild(child, 1), epochIndex, generic);
        }
    }
}

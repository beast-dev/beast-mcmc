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
import dr.evomodel.bigfasttree.BigFastTreeIntervals;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.branchratemodel.TimeVaryingBranchRateModelParser;
import dr.inference.model.AbstractModel;
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

public class TimeVaryingBranchRateModel extends AbstractBranchRateModel
        implements DifferentiableBranchRates, Citable {

    private final Tree tree;
    private final Parameter rates;
    private final EpochTimeProvider epochTimeProvider;

    private boolean nodeRatesKnown;
    private boolean storedNodeRatesKnown;

    private double[] nodeRates;
    private double[] storedNodeRates;

    private final FunctionalForm functionalForm;

    public TimeVaryingBranchRateModel(FunctionalForm.Type type,
                                      Tree tree,
                                      Parameter rates,
                                      Parameter gridPoints) {
        this(type, tree, rates, new EpochTimeProvider.ParameterWrapper(gridPoints));
    }

    public TimeVaryingBranchRateModel(FunctionalForm.Type type,
                                      Tree tree,
                                      Parameter rates,
                                      EpochTimeProvider epochTimeProvider) {

        super(TimeVaryingBranchRateModelParser.PARSER_NAME);

        this.tree = tree;
        this.rates = rates;
        this.epochTimeProvider = epochTimeProvider;

        if (tree instanceof TreeModel) {
            addModel((TreeModel) tree);
        }

        addVariable(rates);
        addModel(epochTimeProvider);

        nodeRates = new double[tree.getNodeCount()];
        storedNodeRates = new double[tree.getNodeCount()];

        functionalForm = type.factory(rates);

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

    public interface EpochTimeProvider extends Model {

        double[] getEpochTimes();

        abstract class AbstractEpochTimeProvider extends AbstractModel implements EpochTimeProvider {

            double[] times;
            boolean timesKnown;

            public AbstractEpochTimeProvider(String name) {
                super(name);
            }

            @Override
            protected void acceptState() { }

            @Override
            public double[] getEpochTimes() {
                if (!timesKnown) {
                    computeTimes();
                    timesKnown = true;
                }
                return times;
            }

            abstract void computeTimes();
        }

        class ParameterWrapper extends AbstractEpochTimeProvider {

            private final Parameter epochTimes;

            public ParameterWrapper(Parameter epochTimes) {
                super("ParameterWrapper");

                this.epochTimes = epochTimes;
                addVariable(epochTimes);
            }

            @Override
            void computeTimes() {
                if (times == null) {
                    times = new double[epochTimes.getDimension() + 1];
                }
                System.arraycopy(epochTimes.getParameterValues(), 0, times, 1,
                        epochTimes.getDimension());
            }

            @Override
            protected void handleModelChangedEvent(Model model, Object object, int index) {
                throw new IllegalArgumentException("Should not be called");
            }

            @Override
            protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
                assert variable == epochTimes;

                timesKnown = false;
                fireModelChanged();
            }

            @Override
            protected void storeState() { }

            @Override
            protected void restoreState() {
                timesKnown = false;
            }
        }

        @SuppressWarnings("unused")
        class IntervalWrapper extends AbstractEpochTimeProvider {

            private final BigFastTreeIntervals intervals;

            public IntervalWrapper(BigFastTreeIntervals intervals) {
                super("IntervalWrapper");

                this.intervals = intervals;
                addModel(intervals);
            }

            @Override
            void computeTimes() {
                int count = intervals.getIntervalCount();

                if (times == null) {
                    times = new double[count];
                }

                for (int i = 0; i < count; ++i) {
                    times[i] = intervals.getIntervalTime(i);
                }
            }

            @Override
            protected void handleModelChangedEvent(Model model, Object object, int index) {
                assert model == intervals;

                timesKnown = false;
                fireModelChanged();
            }

            @Override
            protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
                throw new IllegalArgumentException("Should not be called");
            }

            @Override
            protected void storeState() {
                // TODO
            }

            @Override
            protected void restoreState() {
                // TODO
            }
        }
    }

    public interface FunctionalForm {

        void reset();

        void incrementRate(int epochIndex, double startTime, double endTime);

        double gradientWeight(int epochIndex, double startTime, double endTime, double branchLength);

        double getRateParameter(int epochIndex);

        double rateNumerator();

        enum Type {
            PIECEWISE_CONSTANT("piecewiseConstant") {
                @Override
                FunctionalForm factory(Parameter parameter) {
                    return new PiecewiseConstant(parameter);
                }
            },
            PIECEWISE_LOG_CONSTANT("piecewiseLogConstant") {
                @Override
                FunctionalForm factory(Parameter parameter) {
                    return new PiecewiseLogConstant(parameter);
                }
            };

            private final String name;

            Type(String name) { this.name = name; }

            public String getName() { return name; }

            abstract FunctionalForm factory(Parameter parameter);

            public static Type parse(String string) {
                for (Type type : Type.values()) {
                    if (type.name.equalsIgnoreCase(string)) {
                        return type;
                    }
                }
                throw new IllegalArgumentException("Unknown FunctionalForm.Type");
            }
        }

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
                        getRateParameter(epochIndex);
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

        double[] times = epochTimeProvider.getEpochTimes();
        NodeRef root = tree.getRoot();
        double rootHeight = tree.getNodeHeight(root);

        int epochIndex = times.length - 1;
        while (times[epochIndex] >= rootHeight) {
            --epochIndex;
        }

        traverseTreeByBranchGeneric(times, rootHeight, tree.getChild(root, 0), epochIndex, generic);
        traverseTreeByBranchGeneric(times, rootHeight, tree.getChild(root, 1), epochIndex, generic);
    }

    private void traverseTreeByBranchGeneric(double[] times, double currentHeight, NodeRef child, int epochIndex,
                                             Traversal generic) {

        final double childHeight = tree.getNodeHeight(child);
        final double branchLength = currentHeight - childHeight;
        final int childIndex = getParameterIndexFromNode(child);

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
            traverseTreeByBranchGeneric(times, childHeight, tree.getChild(child, 0), epochIndex, generic);
            traverseTreeByBranchGeneric(times, childHeight, tree.getChild(child, 1), epochIndex, generic);
        }
    }
}

/*
 * AutoCorrelatedBranchRatesDistribution.java
 *
 * Copyright (c) 2002-2019 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
import dr.evomodel.tree.TreeModel;
import dr.inference.distribution.ParametricMultivariateDistributionModel;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.*;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;
import dr.xml.Reportable;

import java.util.Collections;
import java.util.List;

/**
 * @author Marc A. Suchard
 * @author Xiang Ji
 * @author Philippe Lemey
 */
public class AutoCorrelatedBranchRatesDistribution extends AbstractModelLikelihood
        implements GradientWrtParameterProvider, Citable, Reportable {

    private final DifferentiableBranchRates branchRateModel;
    private final ParametricMultivariateDistributionModel distribution;
    private final BranchVarianceScaling scaling;
    private final BranchRateUnits units;

    private final Tree tree;
    private final Parameter rateParameter;

    private boolean incrementsKnown = false;
    private boolean savedIncrementsKnown;

    private boolean likelihoodKnown = false;
    private boolean savedLikelihoodKnown;

    private double logLikelihood;
    private double savedLogLikelihood;

    private double logJacobian;
    private double savedLogJacobian;

    private double priorRateAsIncrement;

    private final int dim;
    private double[] increments;
    private double[] savedIncrements;

    private boolean wrtIncrements;

    public AutoCorrelatedBranchRatesDistribution(String name,
                                                 DifferentiableBranchRates  branchRateModel,
                                                 ParametricMultivariateDistributionModel distribution,
                                                 BranchVarianceScaling scaling,
                                                 boolean takeLogBeforeIncrement, boolean operateOnIncrements) {
        super(name);
        this.branchRateModel = branchRateModel;
        this.distribution = distribution;
        this.scaling = scaling;
        this.units = takeLogBeforeIncrement ? BranchRateUnits.STRICTLY_POSITIVE : BranchRateUnits.REAL_LINE;

        this.tree = branchRateModel.getTree();
        this.rateParameter = branchRateModel.getRateParameter();

        addModel((BranchRateModel) branchRateModel);
        addModel(distribution);

        if (tree instanceof TreeModel) {
            addModel((TreeModel) tree);
        }

        this.dim = branchRateModel.getRateParameter().getDimension();
        this.increments = new double[dim];
        this.savedIncrements = new double[dim];
        this.wrtIncrements = operateOnIncrements;

        if (dim != distribution.getMean().length) {
            throw new RuntimeException("Dimension mismatch in AutoCorrelatedRatesDistribution. " +
                    dim + " != " + distribution.getMean().length);
        }
    }

    public ParametricMultivariateDistributionModel getPrior() {
        return distribution;
    }

    @Override
    public Likelihood getLikelihood() {
        return this;
    }

    @Override
    public Parameter getParameter() {
        return rateParameter;
    }

    @Override
    public int getDimension() {
        return rateParameter.getDimension();
    }

    @Override
    public double[] getGradientLogDensity() {

        double[] gradientWrtIncrements = getGradientWrtIncrements();
        rescaleGradientWrtIncrements(gradientWrtIncrements);

        double[] gradientWrtBranch = new double[dim];
        recurseGradientPreOrder(tree.getRoot(), gradientWrtBranch, gradientWrtIncrements);
        addJacobianTerm(gradientWrtBranch);
        return gradientWrtBranch;
    }

    double[] getGradientWrtIncrements() {

        if (!(distribution instanceof GradientProvider)) {
            throw new RuntimeException("Not yet implemented");
        }

        GradientProvider incrementGradientProvider = (GradientProvider) distribution;
        checkIncrements();
        return incrementGradientProvider.getGradientLogDensity(increments);
    }

    Tree getTree() { return tree; }

    BranchRateUnits getUnits() { return units; }

    BranchVarianceScaling getScaling() { return scaling; }

    DifferentiableBranchRates getBranchRateModel() { return branchRateModel; }

    private void rescaleGradientWrtIncrements(double[] gradientWrtIncrements) {
        for (int i = 0; i < dim; i++) {
            NodeRef node = tree.getNode(i);
            if (!tree.isRoot(node)) {
                    int index = branchRateModel.getParameterIndexFromNode(node);
                gradientWrtIncrements[index] = scaling.rescaleIncrement(
                        gradientWrtIncrements[index], tree.getBranchLength(node));
            }
        }
    }

    private void addJacobianTerm(double[] gradientWrtBranch) {
        for (int i = 0; i < dim; i++) {
            NodeRef node = tree.getNode(i);
            if (!tree.isRoot(node)) {
                int index = branchRateModel.getParameterIndexFromNode(node);
                gradientWrtBranch[index] = units.transformGradient(gradientWrtBranch[index],
                        branchRateModel.getUntransformedBranchRate(tree, node));
            }
        }
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        incrementsKnown = false;
        likelihoodKnown = false;
        fireModelChanged();
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        incrementsKnown = false;
        likelihoodKnown = false;
        fireModelChanged();
    }

    @Override
    protected void storeState() {
        savedIncrementsKnown = incrementsKnown;
        System.arraycopy(increments, 0, savedIncrements, 0, dim);

        savedLikelihoodKnown = likelihoodKnown;
        savedLogLikelihood = logLikelihood;
        savedLogJacobian = logJacobian;
    }

    @Override
    protected void restoreState() {
        incrementsKnown = savedIncrementsKnown;
        double[] tmp = savedIncrements;
        savedIncrements = increments;
        increments = tmp;

        likelihoodKnown = savedLikelihoodKnown;
        logLikelihood = savedLogLikelihood;
        logJacobian = savedLogJacobian;
    }

    @Override
    protected void acceptState() { }

    public double getIncrement(int index) {
        checkIncrements();
        return increments[index];
    }

    @Override
    public Model getModel() {
        return this;
    }

    @Override
    public double getLogLikelihood() {
        if (!likelihoodKnown) {
            logLikelihood = calculateLogLikelihood();
            likelihoodKnown = true;
        }
        return logLikelihood;
    }

    @Override
    public void makeDirty() {
        likelihoodKnown = false;
        incrementsKnown = false;
    }

    @Override
    public Citation.Category getCategory() {
        return null;
    }

    @Override
    public String getDescription() {
        return null;
    }

    @Override
    public List<Citation> getCitations() {
        return Collections.singletonList(CITATION);
    }

    public static Citation CITATION = new Citation(
            new Author[]{
            },
            Citation.Status.IN_PREPARATION
    );

    private void checkIncrements() {
        if (!incrementsKnown) {
            logJacobian =  recursePreOrder(tree.getRoot(), 0); //branchRateModel.getPriorRateAsIncrement(tree));
            incrementsKnown = true;
        }
    }

    private double calculateLogLikelihood() {
        checkIncrements();
        double logLikelihood = distribution.logPdf(increments);
        if (!wrtIncrements) {
            logLikelihood += logJacobian;
        }
        return logLikelihood;
    }

    private double recursePreOrder(NodeRef node, double parentRateAsIncrement) {

        double logJacobian = 0.0;

        if (!tree.isRoot(node)) {
            final double rate = branchRateModel.getUntransformedBranchRate(tree, node);
            final double rateAsIncrement = units.transform(rate);
            final double branchLength = tree.getBranchLength(node);

            logJacobian += units.getTransformLogJacobian(rate) + scaling.getTransformLogJacobian(branchLength);// - branchRateModel.getPriorRateAsIncrement(tree);

            final double rateIncrement = scaling.rescaleIncrement(
                    rateAsIncrement - parentRateAsIncrement, branchLength);

            increments[branchRateModel.getParameterIndexFromNode(node)] = rateIncrement;

            parentRateAsIncrement = rateAsIncrement;
        }

        if (!tree.isExternal(node)) {
            logJacobian += recursePreOrder(tree.getChild(node, 0), parentRateAsIncrement);
            logJacobian += recursePreOrder(tree.getChild(node, 1), parentRateAsIncrement);
        }

        return logJacobian;
    }

    private void recurseGradientPreOrder(NodeRef node,
                                         double[] gradientWrtBranch,
                                         double[] gradientWrtIncrement) {

        int index = branchRateModel.getParameterIndexFromNode(node);

        if (!tree.isRoot(node)) {
            gradientWrtBranch[index] += gradientWrtIncrement[index];
        }

        if (!tree.isExternal(node)) {

            NodeRef child0 = tree.getChild(node, 0);
            NodeRef child1 = tree.getChild(node, 1);

            if (!tree.isRoot(node)) {
                gradientWrtBranch[index] -= gradientWrtIncrement[branchRateModel.getParameterIndexFromNode(child0)];
                gradientWrtBranch[index] -= gradientWrtIncrement[branchRateModel.getParameterIndexFromNode(child1)];
            }

            recurseGradientPreOrder(child0, gradientWrtBranch, gradientWrtIncrement);
            recurseGradientPreOrder(child1, gradientWrtBranch, gradientWrtIncrement);
        }
    }

//    double inverseTransform(double x, NodeRef node) {
//        double branchLength = tree.getBranchLength(node);
//        return units.inverseTransform(scaling.inverseRescaleIncrement(x,branchLength));
//    }

    @Override
    public String getReport() {
        return GradientWrtParameterProvider.getReportAndCheckForError(this,
                0.0, Double.POSITIVE_INFINITY, null);
    }

    public enum BranchRateUnits {

        REAL_LINE("realLine") {
            @Override
            double transform(double x) {
                return x;
            }

            @Override
            double getTransformLogJacobian(double x) {
                return 0.0;
            }

            @Override
            double inverseTransform(double x) {
                return x;
            }

            @Override
            double transformGradient(double gradient, double value) { return gradient; }

            @Override
            double inverseTransformGradient(double gradient, double value) { return gradient; }

            @Override
            boolean needsIncrementCorrection() { return false; }
        },

        STRICTLY_POSITIVE("strictlyPositive") {
            @Override
            double transform(double x) {
                return Math.log(x);
            }

            @Override
            double getTransformLogJacobian(double x) {
                return -Math.log(x);
            }

            @Override
            double inverseTransform(double x) {
                return Math.exp(x);
            }

            @Override
            double transformGradient(double gradient, double value) {
                return (gradient - 1.0) / value;
                // value == r (rate),   \phi = log r
                // gradient_{\phi} == d/d \phi log p(\phi)
                // gradient_{rate} == gradient_{\phi} d \phi / d r + d/d r log Jacobian
                // d \phi / d r == 1/r, d/d r log Jacobian == -1 / r
            }

            @Override
            double inverseTransformGradient(double gradient, double value) {
                return gradient * value;
            }

            @Override
            boolean needsIncrementCorrection() { return true; }
        };

        BranchRateUnits(String name) { this.name = name; }

        public String getName() { return name; }

        private final String name;

        abstract double transform(double x);

        abstract double transformGradient(double gradient, double value);

        abstract double  getTransformLogJacobian(double x);

        abstract double inverseTransform(double x);

        abstract double inverseTransformGradient(double gradient, double value);

        abstract boolean needsIncrementCorrection();
    }

    public enum BranchVarianceScaling {

        NONE("none") {
            @Override
            double rescaleIncrement(double increment, double branchLength) {
                return increment;
            }

            @Override
            double getTransformLogJacobian(double branchLength) { return 0.0; }

            @Override
            double inverseRescaleIncrement(double increment, double branchLength) { return increment; }
        },

        BY_TIME("byTime") {
            @Override
            double rescaleIncrement(double increment, double branchLength) {
                return increment / Math.sqrt(branchLength);
            }

            @Override
            double inverseRescaleIncrement(double increment, double branchLength) {
                return increment * Math.sqrt(branchLength);
            }

            @Override
            double getTransformLogJacobian(double branchLength) {
                return -0.5 * Math.log(branchLength);
            }
        };

        BranchVarianceScaling(String name) {
            this.name = name;
        }

        private final String name;

        abstract double rescaleIncrement(double increment, double branchLength);

        abstract double  inverseRescaleIncrement(double increment, double branchLength);

        abstract double getTransformLogJacobian(double branchLength);

        public String getName() {
            return name;
        }

        public static BranchVarianceScaling parse(String name) {
            for (BranchVarianceScaling scaling : BranchVarianceScaling.values()) {
                if (scaling.getName().equalsIgnoreCase(name)) {
                    return scaling;
                }
            }
            return null;
        }
    }
}

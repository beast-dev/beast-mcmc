/*
 * DiscreteTraitBranchRateGradient.java
 *
 * Copyright (c) 2002-2020 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.treedatalikelihood.discrete;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evolution.tree.TreeTraitProvider;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.DifferentiableBranchRates;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.ProcessSimulation;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.preorder.ProcessSimulationDelegate;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.hmc.HessianWrtParameterProvider;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.math.MultivariateFunction;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;
import dr.xml.Reportable;

import java.util.Collections;
import java.util.List;

import static dr.math.MachineAccuracy.SQRT_EPSILON;


/**
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public class DiscreteTraitBranchRateGradient
        implements GradientWrtParameterProvider, HessianWrtParameterProvider, Reportable, Loggable, Citable {

    protected final TreeDataLikelihood treeDataLikelihood;
    private final TreeTrait treeTraitProvider;
    protected final Tree tree;
    protected final boolean useHessian;
    protected final Parameter rateParameter;
    protected final DifferentiableBranchRates branchRateModel;

    // TODO Refactor / remove code duplication with BranchRateGradient
    // TODO Maybe use:  AbstractBranchRateGradient, DiscreteTraitBranchRateGradient, ContinuousTraitBranchRateGradient
    public DiscreteTraitBranchRateGradient(String traitName,
                                           TreeDataLikelihood treeDataLikelihood,
                                           BeagleDataLikelihoodDelegate likelihoodDelegate,
                                           Parameter rateParameter,
                                           boolean useHessian) {

        assert (treeDataLikelihood != null);

        this.treeDataLikelihood = treeDataLikelihood;
        this.tree = treeDataLikelihood.getTree();
        this.rateParameter = rateParameter;
        this.useHessian = useHessian;

        BranchRateModel brm = treeDataLikelihood.getBranchRateModel();
        this.branchRateModel = (brm instanceof DifferentiableBranchRates) ? (DifferentiableBranchRates) brm : null;

        String name = DiscreteTraitBranchRateDelegate.getName(traitName);
        TreeTrait test = treeDataLikelihood.getTreeTrait(name);

        if (test == null) {
            ProcessSimulationDelegate gradientDelegate = new DiscreteTraitBranchRateDelegate(traitName,
                    treeDataLikelihood.getTree(),
                    likelihoodDelegate);
            TreeTraitProvider traitProvider = new ProcessSimulation(treeDataLikelihood, gradientDelegate);
            treeDataLikelihood.addTraits(traitProvider.getTreeTraits());
        }

        treeTraitProvider = treeDataLikelihood.getTreeTrait(name);
        assert (treeTraitProvider != null);

        int nTraits = treeDataLikelihood.getDataLikelihoodDelegate().getTraitCount();
        if (nTraits != 1) {
            throw new RuntimeException("Not yet implemented for >1 traits");
        }
//        dim = treeDataLikelihood.getDataLikelihoodDelegate().getTraitDim();
    }

    @Override
    public Likelihood getLikelihood() {
        return treeDataLikelihood;
    }

    @Override
    public Parameter getParameter() {
        return rateParameter;
    }

    @Override
    public int getDimension() {
        return getParameter().getDimension();
    }

    public double[] getDiagonalHessianLogDensity() {

        double[] result = new double[tree.getNodeCount() - 1];

        //Do single call to traitProvider with node == null (get full tree)
        double[] diagonalHessian = (double[]) treeDataLikelihood.getTreeTrait(DiscreteTraitBranchRateDelegate.HESSIAN_TRAIT_NAME).getTrait(tree, null);
        double[] gradient = (double[]) treeTraitProvider.getTrait(tree, null);

        int v = 0;
        for (int i = 0; i < tree.getNodeCount(); ++i) {
            final NodeRef node = tree.getNode(i);
            if (!tree.isRoot(node)) {
                final int destinationIndex = getParameterIndexFromNode(node);
                final double chainGradient = getChainGradient(tree, node);
                final double chainSecondDerivative = getChainSecondDerivative(tree, node);
                result[destinationIndex] = diagonalHessian[v] * chainGradient * chainGradient + gradient[v] * chainSecondDerivative;
                v++;
            }
        }

        return result;
    }

    @Override
    public double[][] getHessianLogDensity() {
        throw new RuntimeException("Not yet implemented");
    }

    public double[] getGradientLogDensity() {

        long startTime;
        if (COUNT_TOTAL_OPERATIONS) {
            startTime = System.nanoTime();
        }

        double[] result = new double[tree.getNodeCount() - 1];

        //Do single call to traitProvider with node == null (get full tree)
        double[] gradient = (double[]) treeTraitProvider.getTrait(tree, null);

        int v = 0;
        for (int i = 0; i < tree.getNodeCount(); ++i) {
            final NodeRef node = tree.getNode(i);
            if (!tree.isRoot(node)) {
                final int destinationIndex = getParameterIndexFromNode(node);
                final double nodeResult = gradient[v] * getChainGradient(tree, node);
                result[destinationIndex] = nodeResult;
                v++;
            }
        }

        // TODO Ideally move all chain-ruling into branchRateModel (except branchLengths?)
        result = branchRateModel.updateGradientLogDensity(result, null, 0, gradient.length);

        if (COUNT_TOTAL_OPERATIONS) {
            ++getGradientLogDensityCount;
            long endTime = System.nanoTime();
            totalGradientTime += (endTime - startTime) / 1000000;
        }

        return result;
    }

    protected double getChainGradient(Tree tree, NodeRef node) {
        return tree.getBranchLength(node);
    }

    protected double getChainSecondDerivative(Tree tree, NodeRef node) {
        return 0.0;
    }

    protected int getParameterIndexFromNode(NodeRef node) {
        return (branchRateModel == null) ? node.getNumber() : branchRateModel.getParameterIndexFromNode(node);
    }

    MultivariateFunction numeric1 = new MultivariateFunction() {
        @Override
        public double evaluate(double[] argument) {

            for (int i = 0; i < argument.length; ++i) {
                rateParameter.setParameterValue(i, argument[i]);
            }
            
            return treeDataLikelihood.getLogLikelihood();
        }

        @Override
        public int getNumArguments() {
            return rateParameter.getDimension();
        }

        @Override
        public double getLowerBound(int n) {
            return 0;
        }

        @Override
        public double getUpperBound(int n) {
            return Double.POSITIVE_INFINITY;
        }
    };

    @SuppressWarnings("unused")
    protected boolean valuesAreSufficientlyLarge(double[] vector) {
        for (double x : vector) {
            if (Math.abs(x) < SQRT_EPSILON * 1.2) {
                return false;
            }
        }

        return true;
    }

    @Override
    public String getReport() {

        StringBuilder sb = new StringBuilder();
        if (COUNT_TOTAL_OPERATIONS) {
            sb.append("\n\tgetGradientLogDensityCount = ").append(getGradientLogDensityCount);
            sb.append("\n\taverageGradientTime = ");
            if (getGradientLogDensityCount > 0) {
                sb.append(totalGradientTime / getGradientLogDensityCount);
            } else {
                sb.append("NA");
            }
            sb.append("\n");
        }

        if (CHECK_GRADIENT_IN_REPORT) {
            String message = GradientWrtParameterProvider.getReportAndCheckForError(this, 0.0, Double.POSITIVE_INFINITY, null);

            if (useHessian) {
                message += HessianWrtParameterProvider.getReportAndCheckForError(this, null);
            }

            sb.append(message);
        }

        return  sb.toString();
    }

    private static final boolean CHECK_GRADIENT_IN_REPORT = true;
    protected static final boolean COUNT_TOTAL_OPERATIONS = true;
    private long getGradientLogDensityCount = 0;
    private long totalGradientTime = 0;

    @Override
    public LogColumn[] getColumns() {
        return Loggable.getColumnsFromReport(this, "gradient report");
    }

    @Override
    public Citation.Category getCategory() {
        return Citation.Category.FRAMEWORK;
    }

    @Override
    public String getDescription() {
        return "Using linear-time branch-specific parameter differential calculations";
    }

    @Override
    public List<Citation> getCitations() {
        return Collections.singletonList(CITATION);
    }

    private static final Citation CITATION = new Citation(
            new Author[]{
                    new Author("X", "Ji"),
                    new Author("Z", "Zhang"),
                    new Author("A", "Holbrook"),
                    new Author("A", "Nishimura"),
                    new Author("G", "Beale"),
                    new Author("A", "Rambaut"),
                    new Author("P", "Lemey"),
                    new Author("MA", "Suchard"),
            },
            "Gradients do grow on trees: a linear-time O(N)-dimensional gradient for statistical phylogenetics",
            "Molecular Biology and Evolution",
            Citation.Status.IN_SUBMISSION);
}

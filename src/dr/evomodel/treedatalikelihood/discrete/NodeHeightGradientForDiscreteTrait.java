/*
 * NodeHeightGradientForDiscreteTrait.java
 *
 * Copyright (c) 2002-2017 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
import dr.evolution.tree.TreeTrait;
import dr.evolution.tree.TreeTraitProvider;
import dr.evomodel.branchratemodel.ArbitraryBranchRates;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.ProcessSimulation;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.preorder.AbstractDiscreteTraitDelegate;
import dr.evomodel.treedatalikelihood.preorder.ProcessSimulationDelegate;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.math.MultivariateFunction;
import dr.math.NumericalDerivative;
import dr.xml.Reportable;

import java.util.Arrays;

import static dr.math.MachineAccuracy.SQRT_EPSILON;

/**
 * @author Marc A. Suchard
 * @author Xiang Ji
 */

public class NodeHeightGradientForDiscreteTrait implements GradientWrtParameterProvider, Reportable, Loggable {

    private final TreeDataLikelihood treeDataLikelihood;
    private final TreeTrait treeTraitProvider;
    private final TreeModel tree;
    private final Parameter rateParameter;
    private final ArbitraryBranchRates branchRateModel;
    private Double[] nodeHeights;

    // TODO Refactor / remove code duplication with BranchRateGradient
    // TODO Maybe use:  AbstractBranchRateGradient, DiscreteTraitBranchRateGradient, ContinuousTraitBranchRateGradien

    public NodeHeightGradientForDiscreteTrait(String traitName,
                                              TreeDataLikelihood treeDataLikelihood,
                                              BeagleDataLikelihoodDelegate likelihoodDelegate,
                                              Parameter rateParameter) {
        assert (treeDataLikelihood != null);

        this.treeDataLikelihood = treeDataLikelihood;
        this.tree = (TreeModel) treeDataLikelihood.getTree();
        this.rateParameter = rateParameter;
        this.nodeHeights = new Double[this.tree.getInternalNodeCount()];
        getNodeHeights();

        BranchRateModel brm = treeDataLikelihood.getBranchRateModel();
        this.branchRateModel = (brm instanceof ArbitraryBranchRates) ? (ArbitraryBranchRates) brm : null;
        String name = AbstractDiscreteTraitDelegate.getName(traitName);
        TreeTrait test = treeDataLikelihood.getTreeTrait(name);

        if (test == null) {
            ProcessSimulationDelegate gradientDelegate = new AbstractDiscreteTraitDelegate(traitName,
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

    @Override
    public double[] getGradientLogDensity() {

        double[] result = new double[tree.getInternalNodeCount()];
        Arrays.fill(result, 0.0);

        //Do single call to traitProvider with node == null (get full tree)
        double[] gradient =  (double[]) treeTraitProvider.getTrait(tree, null);

//        if (DEBUG) {
//            System.err.println(new WrappedVector.Raw(gradient));
//        }

        for (int i = 0; i < tree.getInternalNodeCount(); ++i) {
            final NodeRef internalNode = tree.getInternalNode(i);
            if (!tree.isRoot(internalNode)) {
                final int internalNodeNumber = internalNode.getNumber();
                result[i] -= gradient[internalNodeNumber] * branchRateModel.getBranchRate(tree, internalNode);
            }
            for (int j = 0; j < tree.getChildCount(internalNode); ++j){
                NodeRef childNode = tree.getChild(internalNode, j);
                final int childNodeNumber = childNode.getNumber();
                result[i] += gradient[childNodeNumber] * branchRateModel.getBranchRate(tree, childNode);
            }
        }
        return result;
    }

    private double[] getNodeHeights() {
        double[] nodeHeights = new double[tree.getInternalNodeCount()];
        for (int i = 0; i < tree.getInternalNodeCount(); ++i){
            NodeRef internalNode = tree.getInternalNode(i);
            nodeHeights[i] = tree.getNodeHeight(internalNode);
        }
        return nodeHeights;
    }

    private MultivariateFunction numeric1 = new MultivariateFunction() {
        @Override
        public double evaluate(double[] argument) {

            for (int i = 0; i < argument.length; ++i) {
                NodeRef internalNode = tree.getInternalNode(i);
                tree.setNodeHeight(internalNode, argument[i]);
            }

            treeDataLikelihood.makeDirty();
            return treeDataLikelihood.getLogLikelihood();
        }

        @Override
        public int getNumArguments() {
            return tree.getInternalNodeCount();
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

    private boolean valuesAreSufficientlyLarge(double[] vector) {
        for (double x : vector) {
            if (x < SQRT_EPSILON * 1.2) {
                return false;
            }
        }

        return true;
    }


    @Override
    public String getReport() {

        treeDataLikelihood.makeDirty();

        double[] savedValues = getNodeHeights();
        double[] testGradient = null;

        boolean largeEnoughValues = valuesAreSufficientlyLarge(getNodeHeights());

        if (DEBUG && largeEnoughValues) {
            testGradient = NumericalDerivative.gradient(numeric1, getNodeHeights());
        }


        for (int i = 0; i < savedValues.length; ++i) {
            NodeRef internalNode = tree.getInternalNode(i);
            tree.setNodeHeight(internalNode, savedValues[i]);
        }

        StringBuilder sb = new StringBuilder();
        sb.append("Peeling: ").append(new dr.math.matrixAlgebra.Vector(getGradientLogDensity()));
        sb.append("\n");

        if (testGradient != null && largeEnoughValues) {
            sb.append("numeric: ").append(new dr.math.matrixAlgebra.Vector(testGradient));
        } else {
            sb.append("mumeric: too close to 0");
        }
        sb.append("\n");

        treeDataLikelihood.makeDirty();

        return sb.toString();
    }

    private static final boolean DEBUG = true;

    @Override
    public LogColumn[] getColumns() {

        LogColumn[] columns = new LogColumn[1];
        columns[0] = new LogColumn.Default("gradient report", new Object() {
            @Override
            public String toString() {
                return "\n" + getReport();
            }
        });

        return columns;
    }

}

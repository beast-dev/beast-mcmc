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
import dr.evomodel.tree.TreeModel;
import dr.evomodel.tree.TreeParameterModel;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.loggers.Loggable;
import dr.inference.model.Parameter;
import dr.math.MultivariateFunction;
import dr.math.NumericalDerivative;
import dr.xml.Reportable;

import java.util.Arrays;

/**
 * @author Marc A. Suchard
 * @author Xiang Ji
 */

public class NodeHeightGradientForDiscreteTrait extends DiscreteTraitBranchRateGradient
        implements GradientWrtParameterProvider, Reportable, Loggable {

    private final double[] nodeHeights;
    private final TreeModel treeModel;
    protected TreeParameterModel indexHelper;
    private final NodeHeightProxyParameter nodeHeightProxyParameter;


    public NodeHeightGradientForDiscreteTrait(String traitName,
                                              TreeDataLikelihood treeDataLikelihood,
                                              BeagleDataLikelihoodDelegate likelihoodDelegate,
                                              Parameter rateParameter) {
        super(traitName, treeDataLikelihood, likelihoodDelegate, rateParameter, false);

        if (!(treeDataLikelihood.getTree() instanceof TreeModel)) {
            throw new IllegalArgumentException("Must provide a TreeModel");
        }
        this.treeModel = (TreeModel) treeDataLikelihood.getTree();

        this.nodeHeights = new double[tree.getInternalNodeCount()];

        indexHelper = new TreeParameterModel(treeModel, new Parameter.Default(tree.getNodeCount() - 1), false);

        this.nodeHeightProxyParameter = new NodeHeightProxyParameter("internalNodeHeights", treeModel, true);
    }

    @Override
    public Parameter getParameter() {
        return nodeHeightProxyParameter;
    }

    @Override
    public double[] getGradientLogDensity() {

        double[] result = new double[tree.getInternalNodeCount()];
        Arrays.fill(result, 0.0);

        //Do single call to traitProvider with node == null (get full tree)
        double[] gradient = (double[]) treeTraitProvider.getTrait(tree, null);
        

        for (int i = 0; i < tree.getInternalNodeCount(); ++i) {
            final NodeRef internalNode = tree.getInternalNode(i);
            if (!tree.isRoot(internalNode)) {
                final int internalNodeNumber = indexHelper.getParameterIndexFromNodeNumber(internalNode.getNumber());
                result[i] -= gradient[internalNodeNumber] * branchRateModel.getBranchRate(tree, internalNode);
            }
            for (int j = 0; j < tree.getChildCount(internalNode); ++j){
                NodeRef childNode = tree.getChild(internalNode, j);
                final int childNodeNumber = indexHelper.getParameterIndexFromNodeNumber(childNode.getNumber());
                result[i] += gradient[childNodeNumber] * branchRateModel.getBranchRate(tree, childNode);
            }
        }
        return result;
    }

    private double[] getNodeHeights() {
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
                treeModel.setNodeHeight(internalNode, argument[i]);
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
            treeModel.setNodeHeight(internalNode, savedValues[i]);
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


}

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
import dr.evolution.tree.Tree;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.tree.TreeParameterModel;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.loggers.Loggable;
import dr.inference.model.Parameter;
import dr.xml.Reportable;

import java.util.Arrays;

/**
 * @author Marc A. Suchard
 * @author Xiang Ji
 */

public class NodeHeightGradientForDiscreteTrait extends DiscreteTraitBranchRateGradient
        implements GradientWrtParameterProvider, Reportable, Loggable {

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

        indexHelper = new TreeParameterModel(treeModel, new Parameter.Default(tree.getNodeCount() - 1), false);

        this.nodeHeightProxyParameter = new NodeHeightProxyParameter("internalNodeHeights", treeModel, true);
    }

    @Override
    public Parameter getParameter() {
        return nodeHeightProxyParameter;
    }

    protected double getChainGradient(Tree tree, NodeRef node) {
        return branchRateModel.getBranchRate(tree, node);
    }

    @Override
    public double[] getGradientLogDensity() {

        double[] result = new double[tree.getInternalNodeCount()];
        Arrays.fill(result, 0.0);

        double[] gradient = super.getGradientLogDensity();

        for (int i = 0; i < tree.getInternalNodeCount(); ++i) {

            final  NodeRef node = tree.getNode(i + tree.getExternalNodeCount());

            for (int j = 0; j < tree.getChildCount(node); j++) {
                NodeRef childNode = tree.getChild(node, j);
                final int childNodeIndex = indexHelper.getParameterIndexFromNodeNumber(childNode.getNumber());
                result[i] += gradient[childNodeIndex];
            }
            if (!tree.isRoot(node)) {
                final int nodeIndex = indexHelper.getParameterIndexFromNodeNumber(node.getNumber());
                result[i] -= gradient[nodeIndex];
            }
        }
        return result;
    }

    @Override
    public String getReport() {

        String message = GradientWrtParameterProvider.getReportAndCheckForError(this, 0.0, Double.POSITIVE_INFINITY,
                tolerance, smallValueThreshold);

        return message;
    }

    private final double tolerance = 1E-2;
    private final double smallValueThreshold = 1E-3;

    private static final boolean DEBUG = true;


}

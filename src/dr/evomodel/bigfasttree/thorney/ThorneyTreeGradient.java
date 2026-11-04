/*
 * ThorneyTreeGradient.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.evomodel.bigfasttree.thorney;

import dr.evolution.tree.NodeRef;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.tree.TreeParameterModel;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.discrete.NodeHeightProxyParameter;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.xml.Reportable;

/**
 * @author Xiang Ji
 * @author JT McCrone
 */
public class ThorneyTreeGradient implements GradientWrtParameterProvider, Reportable {

    private TreeDataLikelihood likelihood;
    private final NodeHeightProxyParameter nodeHeightProxyParameter;
    private final TreeModel tree;
    private final TreeParameterModel indexHelper;
    private final double[] branchGradient;
    private final ThorneyDataLikelihoodDelegate dataLikelihoodDelegate;
    private final BranchRateModel branchRateModel;
    public ThorneyTreeGradient(TreeDataLikelihood likelihood) {
        this.likelihood = likelihood;
        this.tree = (TreeModel)likelihood.getTree();
        this.dataLikelihoodDelegate = (ThorneyDataLikelihoodDelegate) likelihood.getDataLikelihoodDelegate();
        this.nodeHeightProxyParameter = new NodeHeightProxyParameter("ThorneyTreeGradient.NodeHeightProxyParameter", this.tree, true);
        this.branchGradient = new double[tree.getNodeCount() - 1];
        this.indexHelper = new TreeParameterModel(tree, new Parameter.Default(branchGradient), false);
        this.branchRateModel = likelihood.getBranchRateModel();


    }

    @Override
    public Likelihood getLikelihood() {
        return likelihood;
    }

    @Override
    public Parameter getParameter() {
        return nodeHeightProxyParameter;
    }

    @Override
    public int getDimension() {
        return nodeHeightProxyParameter.getDimension();
    }

    private void calculateBranchGradient() {
        for (int i = 0; i < tree.getNodeCount() - 1; i++) {
            NodeRef node = tree.getNode(indexHelper.getNodeNumberFromParameterIndex(i));
            double time = tree.getBranchLength(node);
            MutationList mutations = dataLikelihoodDelegate.getMutationMap().getMutations(node);
            double rate = branchRateModel.getBranchRate(tree,node);
            branchGradient[i] = dataLikelihoodDelegate.getBranchLengthLikelihoodDelegate().getGradientWrtTime(mutations, time, rate);
        }
    }

    @Override
    public double[] getGradientLogDensity() {
        calculateBranchGradient();
        double[] heightGradient = new double[tree.getInternalNodeCount()];
        for (int i = 0; i < tree.getInternalNodeCount(); ++i) {

            final  NodeRef node = tree.getNode(i + tree.getExternalNodeCount());

            for (int j = 0; j < tree.getChildCount(node); j++) {
                NodeRef childNode = tree.getChild(node, j);
                final int childNodeIndex = indexHelper.getParameterIndexFromNodeNumber(childNode.getNumber());
                heightGradient[i] += branchGradient[childNodeIndex];
            }
            if (!tree.isRoot(node)) {
                final int nodeIndex = indexHelper.getParameterIndexFromNodeNumber(node.getNumber());
                heightGradient[i] -= branchGradient[nodeIndex];
            }
        }
        return heightGradient;
    }


    @Override
    public String getReport() {
        return GradientWrtParameterProvider.getReportAndCheckForError(this, 0, Double.POSITIVE_INFINITY, tolerance);
    }

    private final double tolerance = 1E-3;
}

/*
 * NodeHeightTransform.java
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

import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.coalescent.OldGMRFSkyrideLikelihood;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.continuous.hmc.NodeHeightTransformParser;
import dr.inference.model.Parameter;
import dr.util.Transform;
import dr.xml.Reportable;

/**
 * @author Marc A. Suchard
 * @author Xiang Ji
 */
public class NodeHeightTransform extends Transform.MultivariateTransform implements Reportable {

    private AbstractNodeHeightTransformDelegate nodeHeightTransformDelegate;
    private TreeModel tree;

    public NodeHeightTransform(Parameter nodeHeights,
                               Parameter ratios,
                               TreeModel tree,
                               BranchRateModel branchrateModel) {
        super(nodeHeights.getDimension());
        this.tree = tree;
        if (nodeHeights.getDimension() == tree.getInternalNodeCount()) {
            this.nodeHeightTransformDelegate = new NodeHeightToRatiosFullTransformDelegate(tree, nodeHeights, ratios, branchrateModel);
        } else if (nodeHeights.getDimension() == tree.getInternalNodeCount() - 1) {
            this.nodeHeightTransformDelegate = new NodeHeightToRatiosTransformDelegate(tree, nodeHeights, ratios, branchrateModel);
        } else {
            throw new RuntimeException("Check internal nodeHeight parameter dimentions.");
        }
    }

    public NodeHeightTransform(Parameter nodeHeights,
                               TreeModel tree,
                               OldGMRFSkyrideLikelihood skyrideLikelihood) {
        super(nodeHeights.getDimension());
        this.tree = tree;
        this.nodeHeightTransformDelegate = new NodeHeightToCoalescentIntervalsDelegate(tree, nodeHeights, skyrideLikelihood);
    }

    public Parameter getNodeHeights() {
        return nodeHeightTransformDelegate.getNodeHeights();
    }

    public Parameter getParameter() {
        return nodeHeightTransformDelegate.getParameter();
    }

    @Override
    public double transform(double value) {
        throw new RuntimeException("Should not be called.");
    }

    @Override
    protected double[] transform(double[] values) {
        return nodeHeightTransformDelegate.transform(values);
    }

    @Override
    public double inverse(double value) {
        throw new RuntimeException("Should not be called.");
    }

    @Override
    protected double[] inverse(double[] values) {
        return nodeHeightTransformDelegate.inverse(values);
    }

    @Override
    public boolean isInInteriorDomain(double[] values) {
        return true; //TODO Is it necessary to check entry values ?
    }

    @Override
    public double[] inverse(double[] values, int from, int to, double sum) {
        throw new RuntimeException("Not yet implemented!");
    }

    @Override
    public double[] gradient(double[] values, int from, int to) {
        throw new RuntimeException("Not yet implemented!");
    }

    @Override
    public double[] gradientInverse(double[] values, int from, int to) {
        throw new RuntimeException("Not yet implemented!");
    }

    @Override
    public String getTransformName() {
        return NodeHeightTransformParser.NAME;
    }

    @Override
    protected double getLogJacobian(double[] values) {
        return nodeHeightTransformDelegate.getLogJacobian(values);
    }

    @Override
    protected double[] updateGradientLogDensity(double[] gradient, double[] value) {
        return nodeHeightTransformDelegate.updateGradientLogDensity(gradient, value);
    }

    @Override
    public double[] updateGradientUnWeightedLogDensity(double[] gradient, double[] value, int from, int to) {
        return nodeHeightTransformDelegate.updateGradientUnWeightedLogDensity(gradient, value, from, to);
    }

    @Override
    public double[] updateGradientInverseUnWeightedLogDensity(double[] gradient, double[] value, int from, int to) {
        return nodeHeightTransformDelegate.updateGradientUnWeightedLogDensity(gradient, inverse(value, from, to), from, to);
    }

    @Override
    public String getReport() {
        return nodeHeightTransformDelegate.getReport();
    }

    @Override
    protected double[] getGradientLogJacobianInverse(double[] values) {
        throw new RuntimeException("Not yet implemented!");
    }

    @Override
    public double[][] computeJacobianMatrixInverse(double[] values) {
        throw new RuntimeException("Not yet implemented!");
    }
}

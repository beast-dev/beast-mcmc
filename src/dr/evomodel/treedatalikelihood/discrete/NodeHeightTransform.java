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

    private NodeHeightTransformDelegate nodeHeightTransformDelegate;
    private TreeModel tree;

    public NodeHeightTransform(Parameter nodeHeights,
                               Parameter ratios,
                               TreeModel tree,
                               BranchRateModel branchrateModel) {
        this.tree = tree;
        this.nodeHeightTransformDelegate = new NodeHeightTransformDelegate.Ratios(tree, nodeHeights, ratios, branchrateModel);
    }

    @Override
    public double transform(double value) {
        throw new RuntimeException("Should not be called.");
    }

    @Override
    public double[] transform(double[] values, int from, int to) {
        return nodeHeightTransformDelegate.transform(values, from, to);
    }

    @Override
    public double inverse(double value) {
        throw new RuntimeException("Should not be called.");
    }

    @Override
    public double[] inverse(double[] values, int from, int to) {
        return nodeHeightTransformDelegate.inverse(values, from, to);
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
    public double getLogJacobian(double[] values, int from, int to) {
        throw new RuntimeException("Not yet implemented!");
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

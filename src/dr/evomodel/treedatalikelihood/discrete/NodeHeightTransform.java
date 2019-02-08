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


import dr.evolution.tree.Tree;
import dr.inference.model.Parameter;
import dr.util.Transform;
import dr.xml.Reportable;

/**
 * @author Marc A. Suchard
 * @author Xiang Ji
 */
public class NodeHeightTransform implements Transform, Reportable {

    private Parameter ratios;
    private Parameter nodeHeights;
    private Tree tree;

    public NodeHeightTransform(Parameter nodeHeights, Tree tree) {
        this.nodeHeights = nodeHeights;
        this.tree = tree;
        this.ratios = new Parameter.Default(nodeHeights.getDimension(), 0.5);
    }

    @Override
    public double transform(double value) {
        throw new RuntimeException("Should not be called.");
    }

    @Override
    public double[] transform(double[] values, int from, int to) {
        if (values.length != ratios.getDimension()) {
            throw new RuntimeException("NodeHeightTransform dimension mismatch!");
        }
        updateRatios(values);
        return ratios.getParameterValues();
    }

    private void updateRatios(double[] nodeHeights) {
        // TODO: update ratios here.
    }

    @Override
    public double inverse(double value) {
        throw new RuntimeException("Should not be called.");
    }

    @Override
    public double[] inverse(double[] values, int from, int to) {
        if (values.length != ratios.getDimension()) {
            throw new RuntimeException("NodeHeightTransform dimension mismatch!");
        }
        updateNodeHeights(values);
        return nodeHeights.getParameterValues();
    }

    private void updateNodeHeights(double[] ratios) {
        // TODO: update NodeHeights here.
    }

    @Override
    public double[] inverse(double[] values, int from, int to, double sum) {
        throw new RuntimeException("Not yet implemented!");
    }

    @Override
    public double updateGradientLogDensity(double gradient, double value) {
        throw new RuntimeException("Not yet implemented!");
    }

    @Override
    public double[] updateGradientLogDensity(double[] gradient, double[] value, int from, int to) {
        throw new RuntimeException("Not yet implemented!");
    }

    @Override
    public double updateDiagonalHessianLogDensity(double diagonalHessian, double gradient, double value) {
        throw new RuntimeException("Not yet implemented!");
    }

    @Override
    public double[] updateDiagonalHessianLogDensity(double[] diagonalHessian, double[] gradient, double[] value, int from, int to) {
        throw new RuntimeException("Not yet implemented!");
    }

    @Override
    public double[][] updateHessianLogDensity(double[][] hessian, double[][] transformationHessian, double[] gradient, double[] value, int from, int to) {
        throw new RuntimeException("Not yet implemented!");
    }

    @Override
    public double updateOffdiagonalHessianLogDensity(double offdiagonalHessian, double transformationHessian, double gradientI, double gradientJ, double valueI, double valueJ) {
        throw new RuntimeException("Not yet implemented!");
    }

    @Override
    public double updateGradientInverseUnWeightedLogDensity(double gradient, double value) {
        throw new RuntimeException("Not yet implemented!");
    }

    @Override
    public double[] updateGradientInverseUnWeightedLogDensity(double[] gradient, double[] value, int from, int to) {
        throw new RuntimeException("Not yet implemented!");
    }

    @Override
    public double updateGradientUnWeightedLogDensity(double gradient, double value) {
        throw new RuntimeException("Not yet implemented!");
    }

    @Override
    public double[] updateGradientUnWeightedLogDensity(double[] gradient, double[] value, int from, int to) {
        throw new RuntimeException("Not yet implemented!");
    }

    @Override
    public double gradient(double value) {
        throw new RuntimeException("Not yet implemented!");
    }

    @Override
    public double[] gradient(double[] values, int from, int to) {
        throw new RuntimeException("Not yet implemented!");
    }

    @Override
    public double gradientInverse(double value) {
        throw new RuntimeException("Not yet implemented!");
    }

    @Override
    public double[] gradientInverse(double[] values, int from, int to) {
        throw new RuntimeException("Not yet implemented!");
    }

    @Override
    public String getTransformName() {
        throw new RuntimeException("Not yet implemented!");
    }

    @Override
    public double getLogJacobian(double value) {
        throw new RuntimeException("Not yet implemented!");
    }

    @Override
    public double getLogJacobian(double[] values, int from, int to) {
        throw new RuntimeException("Not yet implemented!");
    }

    @Override
    public boolean isMultivariate() {
        return true;
    }

    @Override
    public String getReport() {
        return "Report ratios here.";
    }
}

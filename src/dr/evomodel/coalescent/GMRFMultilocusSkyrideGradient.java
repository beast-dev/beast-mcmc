/*
 * GMRFMultilocusSkyrideGradient.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.coalescent;

import dr.evolution.tree.NodeRef;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treedatalikelihood.discrete.AbstractNodeHeightTransformDelegate;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.hmc.HessianWrtParameterProvider;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.math.MultivariateFunction;
import dr.math.NumericalDerivative;
import dr.xml.Reportable;

public class GMRFMultilocusSkyrideGradient implements GradientWrtParameterProvider, HessianWrtParameterProvider, Reportable {

    private final GMRFMultilocusSkyrideLikelihood likelihood;
    private final TreeModel tree;
    private final int dim;
    private final Parameter parameter;

    public GMRFMultilocusSkyrideGradient(GMRFMultilocusSkyrideLikelihood likelihood,
                                         TreeModel tree) {
        this.likelihood = likelihood;
        this.tree = tree;
        this.dim = tree.getInternalNodeCount();
        this.parameter = new AbstractNodeHeightTransformDelegate.NodeHeightParameter("allInternalNodes", tree, true);
        tree.addVariable(parameter);
    }


    @Override
    public double[] getDiagonalHessianLogDensity() {
        return new double[dim];
    }



    @Override
    public double[][] getHessianLogDensity() {
        throw new RuntimeException("Not yet implemented!"); // matrix of zeros...
    }

    @Override
    public Likelihood getLikelihood() {
        return likelihood;
    }

    @Override
    public Parameter getParameter() {
        return parameter;
    }

    @Override
    public int getDimension() {
        return dim;
    }

    @Override
    public double[] getGradientLogDensity() {
        return likelihood.getGradientWrtNodeHeights();
    }

    private double[] getGradientWrtNodeHeights() {
        return likelihood.getGradientWrtNodeHeights();
    }

    @Override
    public String getReport() {
        double[] savedValues = getParameter().getParameterValues();
        double[] testGradient = null;

        testGradient = NumericalDerivative.gradient(numeric1, getParameter().getParameterValues());


        for (int i = 0; i < savedValues.length; ++i) {
            NodeRef internalNode = tree.getInternalNode(i);
            tree.setNodeHeightQuietly(internalNode, savedValues[i]);
        }
        tree.pushTreeChangedEvent();

        StringBuilder sb = new StringBuilder();
        sb.append("Peeling: ").append(new dr.math.matrixAlgebra.Vector(getGradientLogDensity()));
        sb.append("\n");

        if (testGradient != null) {
            sb.append("numeric: ").append(new dr.math.matrixAlgebra.Vector(testGradient));
        }
        sb.append("\n");

        return sb.toString();
    }

    private MultivariateFunction numeric1 = new MultivariateFunction() {
        @Override
        public double evaluate(double[] argument) {

            for (int i = 0; i < argument.length; ++i) {
                NodeRef internalNode = tree.getInternalNode(i);
                tree.setNodeHeight(internalNode, argument[i]);
            }

            return likelihood.getLogLikelihood();
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

    private final static Double tolerance = 1E-4;

}

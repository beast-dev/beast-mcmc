/*
 * HyperParameterGradient.java
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

import dr.evolution.tree.MutableTreeModel;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evomodel.tree.TreeParameterModel;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.hmc.HessianWrtParameterProvider;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.math.MultivariateFunction;
import dr.math.NumericalDerivative;
import dr.xml.Reportable;

import static dr.math.MachineAccuracy.SQRT_EPSILON;

/**
 * @author Marc A. Suchard
 * @author Xiang Ji
 */
// TODO: Use this to remove code duplication
public abstract class HyperParameterGradient implements GradientWrtParameterProvider, HessianWrtParameterProvider, Reportable, Loggable {

    final private TreeDataLikelihood treeDataLikelihood;
    final private GradientWrtParameterProvider gradientWrtParameterProvider;
    final private Parameter parameter;
    final private Tree tree;
    final private boolean useHessian;
    final protected TreeParameterModel branchParameter;


    public HyperParameterGradient(TreeDataLikelihood treeDataLikelihood,
                                  GradientWrtParameterProvider gradientWrtParameterProvider,
                                  Parameter parameter,
                                  TreeParameterModel branchParameter,
                                  boolean useHessian) {

        this.treeDataLikelihood = treeDataLikelihood;
        this.gradientWrtParameterProvider = gradientWrtParameterProvider;
        this.parameter = parameter;
        this.useHessian = useHessian;
        this.tree = treeDataLikelihood.getTree();
        this.branchParameter = branchParameter;
    }


    @Override
    public Likelihood getLikelihood() {
        return treeDataLikelihood;
    }

    @Override
    public Parameter getParameter() {
        return parameter;
    }

    @Override
    public int getDimension() {
        return parameter.getDimension();
    }

    @Override
    public double[] getGradientLogDensity() {
        double[] nodeGradients = gradientWrtParameterProvider.getGradientLogDensity();
        if (nodeGradients.length != (tree.getNodeCount() - 1)) {
            throw new RuntimeException("Dimension mismatch!");
        }
        double[] result = new double[getDimension()];
        for (int i = 0; i < branchParameter.getParameterSize(); i++) {
            final NodeRef node = tree.getNode(branchParameter.getNodeNumberFromParameterIndex(i));
            double[] hyperChainGradient = getDifferential(tree, node);
            for (int j = 0; j < result.length; j++) {
                result[j] += nodeGradients[i] * hyperChainGradient[j];
            }
        }
        return result;
    }

    @Override
    public double[] getDiagonalHessianLogDensity() {
        // cannot avoid calculating full hessian in this case, use numerical method for now
        // TODO: maybe add Hessian into BEAGLE ?
        return NumericalDerivative.diagonalHessian(numeric1, parameter.getParameterValues());
    }

    abstract double[] getDifferential(Tree tree, NodeRef node);

    protected MultivariateFunction numeric1 = new MultivariateFunction() {
        @Override
        public double evaluate(double[] argument) {

            for (int i = 0; i < argument.length; ++i) {
                parameter.setParameterValue(i, argument[i]);
            }

//            treeDataLikelihood.makeDirty();
            return treeDataLikelihood.getLogLikelihood();
        }

        @Override
        public int getNumArguments() {
            return parameter.getDimension();
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

//        treeDataLikelihood.makeDirty();

        double[] savedValues = parameter.getParameterValues();
        double[] testGradient = null;
        double[] testHessian = null;

        boolean largeEnoughValues = valuesAreSufficientlyLarge(parameter.getParameterValues());

        if (DEBUG && largeEnoughValues) {
            testGradient = NumericalDerivative.gradient(numeric1, parameter.getParameterValues());
        }

        if (DEBUG && useHessian && largeEnoughValues) {
            testHessian = NumericalDerivative.diagonalHessian(numeric1, parameter.getParameterValues());
        }


        for (int i = 0; i < savedValues.length; ++i) {
            parameter.setParameterValue(i, savedValues[i]);
        }

        StringBuilder sb = new StringBuilder();
        sb.append("Gradient Peeling: ").append(new dr.math.matrixAlgebra.Vector(getGradientLogDensity()));
        sb.append("\n");

        if (testGradient != null && largeEnoughValues) {
            sb.append("Gradient numeric: ").append(new dr.math.matrixAlgebra.Vector(testGradient));
        } else {
            sb.append("Gradient mumeric: too close to 0");
        }
        sb.append("\n");

        if (useHessian) {
            if (largeEnoughValues) {
                sb.append("Hessian Peeling: ").append(new dr.math.matrixAlgebra.Vector(getDiagonalHessianLogDensity()));
                sb.append("\n");
            }

            if (testHessian != null && largeEnoughValues) {
                sb.append("Hessian numeric: ").append(new dr.math.matrixAlgebra.Vector(testHessian));
            } else {
                sb.append("Hessian mumeric: too close to 0");
            }
            sb.append("\n");
        }

        if (COUNT_TOTAL_OPERATIONS) {
            sb.append("\n\tgetGradientLogDensityCount = ").append(getGradientLogDensityCount).append("\n");
            sb.append(treeDataLikelihood.getReport());
        }

        return sb.toString();
    }

    private static final boolean DEBUG = true;

    protected static final boolean COUNT_TOTAL_OPERATIONS = true;
    protected long getGradientLogDensityCount = 0;

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

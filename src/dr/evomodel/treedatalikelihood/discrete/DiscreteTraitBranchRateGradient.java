/*
 * DiscreteTraitBranchRateGradient.java
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
import dr.evolution.tree.TreeTrait;
import dr.evolution.tree.TreeTraitProvider;
import dr.evomodel.branchratemodel.ArbitraryBranchRates;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.ProcessSimulation;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.preorder.AbstractDiscreteTraitDelegate;
import dr.evomodel.treedatalikelihood.preorder.ProcessSimulationDelegate;
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
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public class DiscreteTraitBranchRateGradient
        implements GradientWrtParameterProvider, HessianWrtParameterProvider, Reportable, Loggable {

    protected final TreeDataLikelihood treeDataLikelihood;
    protected final TreeTrait treeTraitProvider;
    protected final Tree tree;
    protected final boolean useHessian;

    //    private final int nTraits;
//    private final int dim;
    protected final Parameter rateParameter;
    protected final ArbitraryBranchRates branchRateModel;

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
        double[] diagonalHessian = (double[]) treeDataLikelihood.getTreeTrait("Hessian").getTrait(tree, null);
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

        double[] result = new double[tree.getNodeCount() - 1];

        //Do single call to traitProvider with node == null (get full tree)
        double[] gradient = (double[]) treeTraitProvider.getTrait(tree, null);

        int v = 0;
        for (int i = 0; i < tree.getNodeCount(); ++i) {
            final NodeRef node = tree.getNode(i);
            if (!tree.isRoot(node)) {
                final int destinationIndex = getParameterIndexFromNode(node);
                final double nodeResult = gradient[v] * getChainGradient(tree, node);
//                if (Double.isNaN(nodeResult) && !Double.isInfinite(treeDataLikelihood.getLogLikelihood())) {
//                    System.err.println("Check Gradient calculation please.");
//                }
                result[destinationIndex] = nodeResult;
                v++;
            }
        }

        if (COUNT_TOTAL_OPERATIONS) {
            ++getGradientLogDensityCount;
        }

        return result;
    }

    protected double getChainGradient(Tree tree, NodeRef node) {
        return tree.getBranchLength(node) / branchRateModel.getBranchRate(tree, node);
    }

    protected double getChainSecondDerivative(Tree tree, NodeRef node) {
        return 0.0;
    }

    protected int getParameterIndexFromNode(NodeRef node) {
        return (branchRateModel == null) ? node.getNumber() : branchRateModel.getParameterIndexFromNode(node);
    }

//    private static final boolean DEBUG = true;

    private MultivariateFunction numeric1 = new MultivariateFunction() {
        @Override
        public double evaluate(double[] argument) {

            for (int i = 0; i < argument.length; ++i) {
                rateParameter.setParameterValue(i, argument[i]);
            }

//            treeDataLikelihood.makeDirty();
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

        double[] savedValues = rateParameter.getParameterValues();
        double[] testGradient = null;
        double[] testHessian = null;

        boolean largeEnoughValues = valuesAreSufficientlyLarge(rateParameter.getParameterValues());

        if (DEBUG && largeEnoughValues) {
            testGradient = NumericalDerivative.gradient(numeric1, rateParameter.getParameterValues());
        }

        if (DEBUG && useHessian && largeEnoughValues) {
            testHessian = NumericalDerivative.diagonalHessian(numeric1, rateParameter.getParameterValues());
        }


        for (int i = 0; i < savedValues.length; ++i) {
            rateParameter.setParameterValue(i, savedValues[i]);
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
            sb.append(treeTraitProvider.toString()).append("\n");
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

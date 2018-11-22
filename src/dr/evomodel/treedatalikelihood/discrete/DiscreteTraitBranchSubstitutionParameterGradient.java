/*
 * DiscreteTraitBranchSubstitutionParameterGradient.java
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

import dr.evolution.tree.*;
import dr.evomodel.branchmodel.ArbitrarySubstitutionParameterBranchModel;
import dr.evomodel.substmodel.DifferentialMassProvider;
import dr.evomodel.treedatalikelihood.discrete.DiscreteTraitBranchSubstitutionParameterDelegate.BranchDifferentialMassProvider;
import dr.evomodel.tree.TreeParameterModel;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.ProcessSimulation;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.preorder.ProcessSimulationDelegate;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.hmc.HessianWrtParameterProvider;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.math.MultivariateFunction;
import dr.math.NumericalDerivative;
import dr.xml.Reportable;

import java.util.List;

import static dr.math.MachineAccuracy.SQRT_EPSILON;

/**
 * @author Marc A. Suchard
 * @author Xiang Ji
 */
public class DiscreteTraitBranchSubstitutionParameterGradient
        implements GradientWrtParameterProvider, HessianWrtParameterProvider, Reportable, Loggable {

    protected final TreeDataLikelihood treeDataLikelihood;
    protected final TreeTrait treeTraitProvider;
    protected final Tree tree;
    protected final ArbitrarySubstitutionParameterBranchModel branchModel;
    protected final boolean useHessian;

    protected final Parameter branchSubstitutionParameter;
    private final TreeParameterModel parameterIndexHelper;

    private static final boolean DEBUG = true;

    protected static final boolean COUNT_TOTAL_OPERATIONS = true;
    protected long getGradientLogDensityCount = 0;

    public DiscreteTraitBranchSubstitutionParameterGradient(String traitName,
                                                            TreeDataLikelihood treeDataLikelihood,
                                                            BeagleDataLikelihoodDelegate likelihoodDelegate,
                                                            Parameter branchSubstitutionParameter,
                                                            ArbitrarySubstitutionParameterBranchModel arbitrarySubstitutionParameterBranchModel,
                                                            List<DifferentialMassProvider> differentialMassProviderList,
                                                            boolean useHessian) {
        this.treeDataLikelihood = treeDataLikelihood;
        this.tree = treeDataLikelihood.getTree();
        this.branchSubstitutionParameter = branchSubstitutionParameter;
        this.branchModel = arbitrarySubstitutionParameterBranchModel;
        this.useHessian = useHessian;
        this.parameterIndexHelper = new TreeParameterModel((MutableTreeModel) tree, new Parameter.Default(tree.getNodeCount() - 1), false);

        String name = DiscreteTraitBranchSubstitutionParameterDelegate.getName(traitName);
        TreeTrait test = treeDataLikelihood.getTreeTrait(name);

        if (test == null) {
            if (!(branchSubstitutionParameter instanceof CompoundParameter)) {
                throw new RuntimeException("Only support CompoundParameter for now.");
            }

            BranchDifferentialMassProvider branchDifferentialMassProvider =
                    new BranchDifferentialMassProvider(parameterIndexHelper, differentialMassProviderList);
            ProcessSimulationDelegate gradientDelegate = new DiscreteTraitBranchSubstitutionParameterDelegate(traitName,
                    treeDataLikelihood.getTree(),
                    likelihoodDelegate,
                    treeDataLikelihood.getBranchRateModel(),
//                    branchModel,
//                    (CompoundParameter) branchSubstitutionParameter,
                    branchDifferentialMassProvider);
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
    public double[] getDiagonalHessianLogDensity() {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public double[][] getHessianLogDensity() {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public Likelihood getLikelihood() {
        return treeDataLikelihood;
    }

    @Override
    public Parameter getParameter() {
        return branchSubstitutionParameter;
    }

    @Override
    public int getDimension() {
        return branchSubstitutionParameter.getDimension();
    }

    @Override
    public double[] getGradientLogDensity() {
        double[] result = new double[tree.getNodeCount() - 1];

        //Do single call to traitProvider with node == null (get full tree)
        double[] gradient = (double[]) treeTraitProvider.getTrait(tree, null);

        int v = 0;
        for (int i = 0; i < tree.getNodeCount(); ++i) {
            final NodeRef node = tree.getNode(i);
            if (!tree.isRoot(node)) {
                final int destinationIndex = parameterIndexHelper.getParameterIndexFromNodeNumber(node.getNumber());
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
        return 1.0;
    }

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

    protected MultivariateFunction numeric1 = new MultivariateFunction() {
        @Override
        public double evaluate(double[] argument) {

            for (int i = 0; i < argument.length; ++i) {
                branchSubstitutionParameter.setParameterValue(i, argument[i]);
            }

//            treeDataLikelihood.makeDirty();
            return treeDataLikelihood.getLogLikelihood();
        }

        @Override
        public int getNumArguments() {
            return branchSubstitutionParameter.getDimension();
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
        double[] savedValues = branchSubstitutionParameter.getParameterValues();
        double[] testGradient = null;
        double[] testHessian = null;

        boolean largeEnoughValues = valuesAreSufficientlyLarge(branchSubstitutionParameter.getParameterValues());

        if (DEBUG && largeEnoughValues) {
            testGradient = NumericalDerivative.gradient(numeric1, branchSubstitutionParameter.getParameterValues());
        }

        if (DEBUG && useHessian && largeEnoughValues) {
            testHessian = NumericalDerivative.diagonalHessian(numeric1, branchSubstitutionParameter.getParameterValues());
        }


        for (int i = 0; i < savedValues.length; ++i) {
            branchSubstitutionParameter.setParameterValue(i, savedValues[i]);
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

        return sb.toString();
    }
}

/*
 * TreeTipGradient.java
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
import dr.evomodel.branchratemodel.ArbitraryBranchRates;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.treedatalikelihood.*;
import dr.evomodel.treedatalikelihood.preorder.AbstractDiscreteTraitDelegate;
import dr.evomodel.treedatalikelihood.preorder.ProcessSimulationDelegate;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.math.MultivariateFunction;
import dr.math.NumericalDerivative;
import dr.math.matrixAlgebra.WrappedVector;
import dr.xml.Reportable;

/**
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public class BranchRateGradientForDiscreteTrait implements GradientWrtParameterProvider, Reportable, Loggable {

    private final TreeDataLikelihood treeDataLikelihood;
    private final TreeTrait treeTraitProvider;
    private final Tree tree;

//    private final int nTraits;
//    private final int dim;
    private final Parameter rateParameter;
    private final ArbitraryBranchRates branchRateModel;

    // TODO Refactor / remove code duplication with BranchRateGradient
    // TODO Maybe use:  AbstractBranchRateGradient, DiscreteTraitBranchRateGradient, ContinuousTraitBranchRateGradien

    public BranchRateGradientForDiscreteTrait(String traitName,
                              TreeDataLikelihood treeDataLikelihood,
                              BeagleDataLikelihoodDelegate likelihoodDelegate,
                              Parameter rateParameter) {

        assert(treeDataLikelihood != null);

        this.treeDataLikelihood = treeDataLikelihood;
        this.tree = treeDataLikelihood.getTree();
        this.rateParameter = rateParameter;

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

    @Override
    public double[] getGradientLogDensity() {

        double[] result = new double[rateParameter.getDimension()];

        //Do single call to traitProvider with node == null (get full tree)
        double[] gradient =  (double[]) treeTraitProvider.getTrait(tree, null);

        System.err.println(new WrappedVector.Raw(gradient));

        int v =0;
        for (int i = 0; i < tree.getNodeCount(); ++i) {
            final NodeRef node = tree.getNode(i);
            if (!tree.isRoot(node)) {
                final int destinationIndex = getParameterIndexFromNode(node);
                result[destinationIndex] = gradient[v++];
            }
        }

        return result;
    }

    private int getParameterIndexFromNode(NodeRef node) {
        return (branchRateModel == null) ? node.getNumber() : branchRateModel.getParameterIndexFromNode(node);
    }

//    private static final boolean DEBUG = true;

    private MultivariateFunction numeric1 = new MultivariateFunction() {
        @Override
        public double evaluate(double[] argument) {

            for (int i = 0; i < argument.length; ++i) {
                rateParameter.setParameterValue(i, argument[i]);
            }

            treeDataLikelihood.makeDirty();
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

    @Override
    public String getReport() {
        double[] savedValues = rateParameter.getParameterValues();
        double[] testGradient;

        if (DEBUG) {
            testGradient = NumericalDerivative.gradient(numeric1, rateParameter.getParameterValues());
        }


        for (int i = 0; i < savedValues.length; ++i) {
            rateParameter.setParameterValue(i, savedValues[i]);
        }

        StringBuilder sb = new StringBuilder();
        sb.append("Peeling: ").append(new dr.math.matrixAlgebra.Vector(getGradientLogDensity()));
        sb.append("\n");
        
        if (DEBUG) {
            sb.append("numeric: ").append(new dr.math.matrixAlgebra.Vector(testGradient));
            sb.append("\n");
        }

        return sb.toString();
    }

    private static final boolean DEBUG = false;

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

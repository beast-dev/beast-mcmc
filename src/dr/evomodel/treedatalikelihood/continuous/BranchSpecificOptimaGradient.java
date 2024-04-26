/*
 * BranchSpecificOptimaGradient.java
 *
 * Copyright (c) 2002-2019 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.treedatalikelihood.continuous;


import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evomodel.branchratemodel.ArbitraryBranchRates;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.preorder.BranchConditionalDistributionDelegate;
import dr.evomodel.treedatalikelihood.preorder.BranchSufficientStatistics;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.math.MultivariateFunction;
import dr.math.NumericalDerivative;
import dr.xml.Reportable;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Alexander Fisher
 */
public class BranchSpecificOptimaGradient implements GradientWrtParameterProvider, Reportable {

    private TreeDataLikelihood treeDataLikelihood;
    private ContinuousTraitGradientForBranch branchProvider;
    private final TreeTrait<List<BranchSufficientStatistics>> treeTraitProvider;
    private final int numBranches;
    private final int numTraits;
    private final Tree tree;
    private final int dimension;
    private Parameter parameter;

    public BranchSpecificOptimaGradient(String traitName, TreeDataLikelihood treeDataLikelihood, ContinuousDataLikelihoodDelegate likelihoodDelegate, ContinuousTraitGradientForBranch branchProvider, CompoundParameter optimaParameter) {
        this.branchProvider = branchProvider;
        this.treeDataLikelihood = treeDataLikelihood;
        this.tree = treeDataLikelihood.getTree();
        this.numTraits = optimaParameter.getParameterCount();
        this.numBranches = treeDataLikelihood.getTree().getNodeCount() - 1;
        this.parameter = optimaParameter;
        this.dimension = parameter.getDimension();

        String bcdName = BranchConditionalDistributionDelegate.getName(traitName);
        if (treeDataLikelihood.getTreeTrait(bcdName) == null) {
            likelihoodDelegate.addBranchConditionalDensityTrait(traitName);
        }

        TreeTrait<List<BranchSufficientStatistics>> unchecked = treeDataLikelihood.getTreeTrait(bcdName);
        treeTraitProvider = unchecked;

        assert (treeTraitProvider != null);

        getGradientLogDensity();
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
        return dimension;
    }

    public double[] getGradientLogDensity() {
        double[] gradient = new double[numBranches * numTraits];
        double[] gradientAtBranchI = new double[numTraits];

        int index = 0;

        for (int i = 0; i < numBranches; i++) {
            final NodeRef node = tree.getNode(i);
            List<BranchSufficientStatistics> statisticsForNode = treeTraitProvider.getTrait(tree, node);
            gradientAtBranchI = branchProvider.getGradientForBranch(statisticsForNode.get(0), node);
            for (int traitDim = 0; traitDim < numTraits; traitDim++) {
                gradient[i + (numBranches * traitDim)] = gradientAtBranchI[traitDim];
                index = index + 1;
            }
        }
        return gradient;
    }

    private MultivariateFunction numeric1 = new MultivariateFunction() {
        @Override
        public double evaluate(double[] argument) {

            for (int i = 0; i < dimension; i++) {
                parameter.setParameterValue(i, argument[i]);
            }

            treeDataLikelihood.makeDirty();
            return treeDataLikelihood.getLogLikelihood();
        }

        @Override
        public int getNumArguments() {
            return dimension;
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

    public double[] getNumericalGradient() {
        double[] savedValues = parameter.getParameterValues();

        double[] testGradient = NumericalDerivative.gradient(numeric1, parameter.getParameterValues());

        for (int i = 0; i < dimension; i++) {
            parameter.setParameterValue(i, savedValues[i]);
        }

        return testGradient;
    }

    @Override
    public String getReport() {
        double[] testGradient = getNumericalGradient();

        StringBuilder sb = new StringBuilder();
        sb.append("peeling: ").append(new dr.math.matrixAlgebra.Vector(getGradientLogDensity()));
        sb.append("\n");
        sb.append("numeric: ").append(new dr.math.matrixAlgebra.Vector(testGradient));
        sb.append("\n");

        return sb.toString();
    }
}

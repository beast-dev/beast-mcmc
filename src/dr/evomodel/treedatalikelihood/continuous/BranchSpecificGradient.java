/*
 * BranchSpecificGradient.java
 *
 * Copyright (c) 2002-2018 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.preorder.BranchConditionalDistributionDelegate;
import dr.evomodel.treedatalikelihood.preorder.BranchSufficientStatistics;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.math.MultivariateFunction;
import dr.math.NumericalDerivative;
import dr.xml.Reportable;

import java.util.List;

/**
 * @author Marc A. Suchard
 * @author Paul Bastide
 */
public class BranchSpecificGradient implements GradientWrtParameterProvider, Reportable, Loggable {

    private final TreeDataLikelihood treeDataLikelihood;
    private final TreeTrait<List<BranchSufficientStatistics>> treeTraitProvider;
    private final Tree tree;
    private final int nTraits;
    //    private final int dim;
    private final Parameter parameter;
    //    private final ArbitraryBranchRates branchRateModel;
    private final ContinuousTraitGradientForBranch branchProvider;

    public BranchSpecificGradient(String traitName,
                                  TreeDataLikelihood treeDataLikelihood,
                                  ContinuousDataLikelihoodDelegate likelihoodDelegate,
                                  ContinuousTraitGradientForBranch branchProvider,
                                  Parameter parameter) {

        assert (treeDataLikelihood != null);

        this.treeDataLikelihood = treeDataLikelihood;
        this.tree = treeDataLikelihood.getTree();
        this.parameter = parameter;
        this.branchProvider = branchProvider;

        // TODO Move into different constructor / parser
        String bcdName = BranchConditionalDistributionDelegate.getName(traitName);
        if (treeDataLikelihood.getTreeTrait(bcdName) == null) {
            likelihoodDelegate.addBranchConditionalDensityTrait(traitName);
        }

        @SuppressWarnings("unchecked")
        TreeTrait<List<BranchSufficientStatistics>> unchecked = treeDataLikelihood.getTreeTrait(bcdName);
        treeTraitProvider = unchecked;

        assert (treeTraitProvider != null);

        nTraits = treeDataLikelihood.getDataLikelihoodDelegate().getTraitCount();
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
        return parameter;
    }

    @Override
    public int getDimension() {
        return getParameter().getDimension();
    }

    @Override
    public double[] getGradientLogDensity() {

        int dimGradient = branchProvider.getDimension();

        double[] result = new double[parameter.getDimension()];

        // TODO Do single call to traitProvider with node == null (get full tree)
//        List<BranchSufficientStatistics> statisticsForTree = (List<BranchSufficientStatistics>)
//                treeTraitProvider.getTrait(tree, null);

        for (int i = 0; i < tree.getNodeCount(); ++i) {
            final NodeRef node = tree.getNode(i);

//            if (!tree.isRoot(node)) {

            List<BranchSufficientStatistics> statisticsForNode = treeTraitProvider.getTrait(tree, node);

            assert (statisticsForNode.size() == nTraits);

            double[] gradient;
//                for (int trait = 0; trait < nTraits; ++trait) { // TODO deal with several traits
            gradient = branchProvider.getGradientForBranch(statisticsForNode.get(0), node);
//                }

            final int destinationIndex = getParameterIndexFromNode(node);
            assert (destinationIndex != -1);
            for (int j = 0; j < dimGradient; j++) {
                result[destinationIndex * dimGradient + j] += gradient[j];
            }
//            }
        }

        return result;
    }

    private int getParameterIndexFromNode(NodeRef node) {
        return branchProvider.getParameterIndexFromNode(node);
//        return (branchRateModel == null) ? node.getNumber() : branchRateModel.getParameterIndexFromNode(node);
    }


    private MultivariateFunction numeric1 = new MultivariateFunction() {
        @Override
        public double evaluate(double[] argument) {

            for (int i = 0; i < argument.length; ++i) {
                parameter.setParameterValue(i, argument[i]);
            }

            treeDataLikelihood.makeDirty();
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

    @Override
    public String getReport() {

        double[] savedValues = parameter.getParameterValues();
        double[] testGradient = NumericalDerivative.gradient(numeric1, parameter.getParameterValues());
        for (int i = 0; i < savedValues.length; ++i) {
            parameter.setParameterValue(i, savedValues[i]);
        }

        StringBuilder sb = new StringBuilder();
        sb.append("Peeling: ").append(new dr.math.matrixAlgebra.Vector(getGradientLogDensity()));
        sb.append("\n");
        sb.append("numeric: ").append(new dr.math.matrixAlgebra.Vector(testGradient));
        sb.append("\n");

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


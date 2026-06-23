/*
 * CanonicalOUNodeHeightGradient.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.evomodel.treedatalikelihood.continuous;

import dr.evomodel.branchratemodel.ArbitraryBranchRates;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.DefaultBranchRateModel;
import dr.evomodel.branchratemodel.StrictClockBranchRates;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.discrete.NodeHeightProxyParameter;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.math.MultivariateFunction;
import dr.math.NumericalDerivative;
import dr.xml.Reportable;

/**
 * Internal-node-height gradient for canonical OU tree likelihoods with height-independent branch rates.
 */
public final class CanonicalOUNodeHeightGradient
        implements GradientWrtParameterProvider, Reportable, Loggable {

    private final TreeDataLikelihood treeDataLikelihood;
    private final ContinuousDataLikelihoodDelegate likelihoodDelegate;
    private final TreeModel tree;
    private final BranchRateModel branchRateModel;
    private final NodeHeightProxyParameter nodeHeightProxyParameter;

    public CanonicalOUNodeHeightGradient(final TreeDataLikelihood treeDataLikelihood,
                                         final ContinuousDataLikelihoodDelegate likelihoodDelegate) {
        if (!(treeDataLikelihood.getTree() instanceof TreeModel)) {
            throw new IllegalArgumentException("Canonical OU node-height gradients require a TreeModel.");
        }
        if (!likelihoodDelegate.usesCanonicalOULikelihood()) {
            throw new IllegalArgumentException(
                    "Canonical OU node-height gradient requires canonical OU likelihood mode.");
        }
        this.treeDataLikelihood = treeDataLikelihood;
        this.likelihoodDelegate = likelihoodDelegate;
        this.tree = (TreeModel) treeDataLikelihood.getTree();
        this.branchRateModel = treeDataLikelihood.getBranchRateModel();
        validateBranchRateModel(branchRateModel);
        this.nodeHeightProxyParameter =
                new NodeHeightProxyParameter("internalNodeHeights", tree, true);
    }

    public static boolean supportsBranchRateModel(final BranchRateModel branchRateModel) {
        return branchRateModel == null
                || branchRateModel instanceof DefaultBranchRateModel
                || branchRateModel instanceof StrictClockBranchRates
                || branchRateModel instanceof ArbitraryBranchRates;
    }

    private static void validateBranchRateModel(final BranchRateModel branchRateModel) {
        if (!supportsBranchRateModel(branchRateModel)) {
            throw new IllegalArgumentException(
                    "Canonical OU node-height gradients require branch rates that do not depend on node heights. "
                            + "Unsupported branch-rate model: " + branchRateModel.getClass().getName());
        }
    }

    @Override
    public Likelihood getLikelihood() {
        return treeDataLikelihood;
    }

    @Override
    public Parameter getParameter() {
        return nodeHeightProxyParameter;
    }

    @Override
    public int getDimension() {
        return getParameter().getDimension();
    }

    @Override
    public double[] getGradientLogDensity() {
        treeDataLikelihood.makeDirty();
        final double[] effectiveTimeGradient =
                likelihoodDelegate.getCanonicalOUGradientAdapter()
                        .getEffectiveBranchTimeGradient(tree.getNodeCount());
        final double[] gradient = new double[getDimension()];
        CanonicalOUBranchTimeGradient.fillNodeHeightGradient(
                tree,
                branchRateModel,
                likelihoodDelegate.getRateTransformation(),
                effectiveTimeGradient,
                gradient);
        return gradient;
    }

    public double[] getNumericalGradient() {
        final double[] savedValues = getParameter().getParameterValues();
        final double[] testGradient =
                NumericalDerivative.gradient(numericLikelihood, savedValues.clone());
        for (int i = 0; i < savedValues.length; ++i) {
            getParameter().setParameterValueQuietly(i, savedValues[i]);
        }
        getParameter().fireParameterChangedEvent();
        treeDataLikelihood.makeDirty();
        return testGradient;
    }

    @Override
    public String getReport() {
        StringBuilder sb = new StringBuilder();
        sb.append("Canonical OU node-height gradient: ")
                .append(new dr.math.matrixAlgebra.Vector(getGradientLogDensity()));
        sb.append("\n");
        sb.append("numeric: ")
                .append(new dr.math.matrixAlgebra.Vector(getNumericalGradient()));
        sb.append("\n");
        return sb.toString();
    }

    @Override
    public LogColumn[] getColumns() {
        return Loggable.getColumnsFromReport(this, "CanonicalOUNodeHeightGradient report");
    }

    private final MultivariateFunction numericLikelihood = new MultivariateFunction() {
        @Override
        public double evaluate(final double[] argument) {
            for (int i = 0; i < argument.length; ++i) {
                getParameter().setParameterValueQuietly(i, argument[i]);
            }
            getParameter().fireParameterChangedEvent();
            treeDataLikelihood.makeDirty();
            return treeDataLikelihood.getLogLikelihood();
        }

        @Override
        public int getNumArguments() {
            return getDimension();
        }

        @Override
        public double getLowerBound(final int n) {
            return 0.0;
        }

        @Override
        public double getUpperBound(final int n) {
            return Double.POSITIVE_INFINITY;
        }
    };
}

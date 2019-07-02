/*
 * AutoCorrelatedGradientWrtIncrements.java
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

package dr.evomodel.branchratemodel;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.xml.Reportable;

/**
 * @author Marc A. Suchard
 * @author Alex Fisher
 * @author Xiang Ji
 */
public class BranchRateGradientWrtIncrements implements GradientWrtParameterProvider, Reportable {

    private final GradientWrtParameterProvider rateGradientProvider;
    private final AutoCorrelatedGradientWrtIncrements priorGradientProvider;

    private final ArbitraryBranchRates branchRates;
    private final Tree tree;

    private final AutoCorrelatedBranchRatesDistribution.BranchVarianceScaling scaling;
    private final AutoCorrelatedBranchRatesDistribution.BranchRateUnits units;

    public BranchRateGradientWrtIncrements(GradientWrtParameterProvider rateGradientProvider,
                                           AutoCorrelatedGradientWrtIncrements priorGradientProvider) {
        this.rateGradientProvider = rateGradientProvider;
        this.priorGradientProvider = priorGradientProvider;

        AutoCorrelatedBranchRatesDistribution distribution = priorGradientProvider.getDistribution();

        this.branchRates = distribution.getBranchRateModel();
        this.tree = distribution.getTree();
        this.scaling = distribution.getScaling();
        this.units = distribution.getUnits();
    }

    @Override
    public Likelihood getLikelihood() {
        return rateGradientProvider.getLikelihood();
    }

    @Override
    public Parameter getParameter() {
        return priorGradientProvider.getParameter();
    }

    @Override
    public int getDimension() {
        return priorGradientProvider.getDimension();
    }

    @Override
    public double[] getGradientLogDensity() {

        double[] gradientWrtIncrements = rateGradientProvider.getGradientLogDensity();
        double[] gradientWrtRates = new double[gradientWrtIncrements.length];

        recursePostOrderToAccumulateGradient(tree.getRoot(), gradientWrtRates, gradientWrtIncrements);

        return gradientWrtRates;
    }

    private double recursePostOrderToAccumulateGradient(NodeRef node,
                                                        double[] gradientWrtRates, double[] gradientWrtIncrements) {

        // On STRICTLY_POSITIVE scale, log-likelihood includes log-Jacobian (\sum_{increments} -> rate)

        double gradientForNode = 0.0;

        if (!tree.isExternal(node)) {
            gradientForNode += recursePostOrderToAccumulateGradient(tree.getChild(node, 0),
                    gradientWrtRates, gradientWrtIncrements);
            gradientForNode += recursePostOrderToAccumulateGradient(tree.getChild(node, 1),
                    gradientWrtRates, gradientWrtIncrements);
        }

        if (!tree.isRoot(node)) {
            int index = branchRates.getParameterIndexFromNode(node);
            gradientForNode += units.inverseTransformGradient(
                    gradientWrtIncrements[index], branchRates.getUntransformedBranchRate(tree, node));
            gradientWrtRates[index] = scaling.inverseRescaleIncrement(gradientForNode, tree.getBranchLength(node));
        }

        return gradientForNode;
    }

    @Override
    public String getReport() {
        return GradientWrtParameterProvider.getReportAndCheckForError(this,
                Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY, null);
    }
}

/*
 * BranchRateGradientWrtIncrements.java
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
import dr.inference.model.CompoundParameter;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.xml.Reportable;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Marc A. Suchard
 * @author Alexander Fisher
 * @author Xiang Ji
 */
public class BranchRateGradientWrtIncrements implements GradientWrtParameterProvider, Reportable {

    private final GradientWrtParameterProvider rateGradientProvider;
    private List<AutoCorrelatedGradientWrtIncrements> priorGradientProvider;

    private List<DifferentiableBranchRates> branchRates = new ArrayList<>();
    private final Tree tree;

    private List<AutoCorrelatedBranchRatesDistribution.BranchVarianceScaling> scaling = new ArrayList<>();
    private List<AutoCorrelatedBranchRatesDistribution.BranchRateUnits> units = new ArrayList<>();

    private final int listDim;
    private final int paramDim;

    private final Parameter parameter;

    public BranchRateGradientWrtIncrements(GradientWrtParameterProvider rateGradientProvider,
                                           List<AutoCorrelatedGradientWrtIncrements> priorGradientProvider) {
        this.rateGradientProvider = rateGradientProvider;
        this.priorGradientProvider = priorGradientProvider;

        this.listDim = priorGradientProvider.size();
        List<AutoCorrelatedBranchRatesDistribution> distribution = new ArrayList<>();
        AutoCorrelatedBranchRatesDistribution dist;

        CompoundParameter param = new CompoundParameter(null);

        for (int i = 0; i < listDim; i++) {
            dist = priorGradientProvider.get(i).getDistribution();
            distribution.add(dist);
            branchRates.add(dist.getBranchRateModel());
            scaling.add(dist.getScaling());
            units.add(dist.getUnits());
            param.addParameter(priorGradientProvider.get(i).getParameter());
        }

        this.parameter = param;
        this.paramDim = distribution.get(0).getParameter().getDimension();
//        this.branchRates = distribution.getBranchRateModel();
        this.tree = distribution.get(0).getTree();
//        this.scaling = distribution.getScaling();
//        this.units = distribution.getUnits();
    }

    @Override
    public Likelihood getLikelihood() {
        return rateGradientProvider.getLikelihood();
    }

    @Override
    public Parameter getParameter() {
//        return priorGradientProvider.getParameter();
        return parameter;
    }

    @Override
    public int getDimension() {
//        return priorGradientProvider.getDimension();
        return parameter.getDimension();
    }

    @Override
    public double[] getGradientLogDensity() {

        double[] fullGradientWrtRates = rateGradientProvider.getGradientLogDensity();
        double[] fullGradientWrtIncrements = new double[fullGradientWrtRates.length];

        double[] gradientWrtRates = new double[paramDim]; // rateGradientProvider.getGradientLogDensity();
        double[] gradientWrtIncrements = new double[gradientWrtRates.length];

        for (int paramIndex = 0; paramIndex < listDim; paramIndex++) {
            System.arraycopy(fullGradientWrtRates, paramIndex * paramDim, gradientWrtRates, 0, paramDim);
            recursePostOrderToAccumulateGradient(tree.getRoot(), gradientWrtRates, gradientWrtIncrements, paramIndex);
            System.arraycopy(gradientWrtIncrements, 0, fullGradientWrtIncrements, paramIndex * paramDim, paramDim);
        }

        return fullGradientWrtIncrements;
    }

    private double recursePostOrderToAccumulateGradient(NodeRef node,
                                                        double[] gradientWrtRates, double[] gradientWrtIncrements, int paramIndex) {

        // On STRICTLY_POSITIVE scale, log-likelihood includes log-Jacobian (\sum_{increments} -> rate)

        double gradientForNode = 0.0;
//        for (int i = 0; i < listDim; i++) {
        if (!tree.isExternal(node)) {
            gradientForNode += recursePostOrderToAccumulateGradient(tree.getChild(node, 0),
                    gradientWrtRates, gradientWrtIncrements, paramIndex);
            gradientForNode += recursePostOrderToAccumulateGradient(tree.getChild(node, 1),
                    gradientWrtRates, gradientWrtIncrements, paramIndex);
        }

        if (!tree.isRoot(node)) {
            int index = branchRates.get(paramIndex).getParameterIndexFromNode(node);
            gradientForNode += units.get(paramIndex).inverseTransformGradient(
                    gradientWrtRates[index], branchRates.get(paramIndex).getUntransformedBranchRate(tree, node));
//                        gradientWrtRates[index + (paramIndex * paramDim)], branchRates.get(paramIndex).getUntransformedBranchRate(tree, node));
//                gradientWrtIncrements[index + (paramIndex * paramDim)] = scaling.get(paramIndex).inverseRescaleIncrement(gradientForNode, tree.getBranchLength(node));
            gradientWrtIncrements[index] = scaling.get(paramIndex).inverseRescaleIncrement(gradientForNode, tree.getBranchLength(node));

        }
//        }

        return gradientForNode;
    }

    @Override
    public String getReport() {
        return GradientWrtParameterProvider.getReportAndCheckForError(this,
                Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY, null);
    }
}

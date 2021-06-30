/*
 * JointBayesianBridgeStatistics.java
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

package dr.inference.distribution.shrinkage;

import dr.evomodel.branchratemodel.AutoCorrelatedBranchRatesDistribution;
import dr.inference.model.Parameter;

import java.util.List;

/**
 * @author Alexander Fisher
 */

public class JointBayesianBridgeStatistics implements BayesianBridgeStatisticsProvider {

    private List<AutoCorrelatedBranchRatesDistribution> incProviderList;
    private BayesianBridgeDistributionModel bayesianBridgeDistribution;

    private final int dim;

    private double[] branchSpecificFunction;

    public JointBayesianBridgeStatistics(List<AutoCorrelatedBranchRatesDistribution> incProviderList) {
        this.incProviderList = incProviderList;
        this.bayesianBridgeDistribution = (BayesianBridgeDistributionModel) incProviderList.get(0).getPrior();

        dim = incProviderList.get(0).getDimension();

        this.branchSpecificFunction = new double[dim];
        computeBranchSpecificFunction();
    }

    private void computeBranchSpecificFunction() {
        double totalIncSquared;
        for (int i = 0; i < dim; i++) {
            totalIncSquared = 0;
            for (int j = 0; j < incProviderList.size(); j++) {
                totalIncSquared += Math.pow(incProviderList.get(j).getIncrement(i), 2);
            }
            branchSpecificFunction[i] = Math.sqrt(totalIncSquared);
        }
    }

    @Override
    public double getCoefficient(int i) {
        return branchSpecificFunction[i];
    }

    @Override
    public Parameter getGlobalScale() {
        return bayesianBridgeDistribution.getGlobalScale();
    }

    @Override
    public Parameter getLocalScale() {
        return bayesianBridgeDistribution.getLocalScale();
    }

    @Override
    public Parameter getExponent() {
        return bayesianBridgeDistribution.getExponent();
    }

    @Override
    public int getDimension() {
        return dim;
    }
}

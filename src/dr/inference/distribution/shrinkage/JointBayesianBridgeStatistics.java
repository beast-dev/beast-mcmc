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

    private final List<BayesianBridgeStatisticsProvider> providers;

    public JointBayesianBridgeStatistics(List<BayesianBridgeStatisticsProvider> providers) {

        this.providers = providers;

        BayesianBridgeStatisticsProvider base = providers.get(0);
        for (int i = 1; i < providers.size(); ++i) {
            BayesianBridgeStatisticsProvider next = providers.get(i);
            if (base.getDimension() != next.getDimension() ||
                    base.getExponent() != next.getExponent() ||
                    base.getGlobalScale() != next.getGlobalScale() ||
                    base.getLocalScale() != next.getLocalScale()) {
                throw new IllegalArgumentException("All Bayesian bridges must be the same");
            }
        }
    }

    @Override
    public double getCoefficient(int i) {

        double squaredSum = 0.0;
        for (BayesianBridgeStatisticsProvider p : providers) {
            squaredSum += p.getCoefficient(i) * p.getCoefficient(i);
        }

        return Math.sqrt(squaredSum);
    }

    @Override
    public Parameter getGlobalScale() {
        return providers.get(0).getGlobalScale();
    }

    @Override
    public Parameter getLocalScale() {
        return providers.get(0).getLocalScale();
    }

    @Override
    public Parameter getExponent() {
        return providers.get(0).getExponent();
    }

    @Override
    public int getDimension() {
        return providers.get(0).getDimension();
    }
}

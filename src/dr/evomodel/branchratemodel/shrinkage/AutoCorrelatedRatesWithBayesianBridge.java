/*
 * AutoCorrelatedRatesWithBayesianBridge.java
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

package dr.evomodel.branchratemodel.shrinkage;

import dr.evomodel.branchratemodel.AutoCorrelatedBranchRatesDistribution;
import dr.inference.distribution.shrinkage.BayesianBridgeDistributionModel;
import dr.inference.distribution.shrinkage.BayesianBridgeStatisticsProvider;
import dr.inference.model.Parameter;

/**
 * @author Marc A. Suchard
 */

@Deprecated
public class AutoCorrelatedRatesWithBayesianBridge implements BayesianBridgeStatisticsProvider {

    private final AutoCorrelatedBranchRatesDistribution rates;
    private final BayesianBridgeDistributionModel prior;

    public AutoCorrelatedRatesWithBayesianBridge(AutoCorrelatedBranchRatesDistribution rates,
                                                 BayesianBridgeDistributionModel prior) {
        this.rates = rates;
        this.prior = prior;
    }

    @Override
    public double getCoefficient(int i) {
        return rates.getIncrement(i);
    }

    @Override
    public Parameter getGlobalScale() {
        return prior.getGlobalScale();
    }

    @Override
    public Parameter getLocalScale() {
        return prior.getLocalScale();
    }

    @Override
    public Parameter getExponent() {
        return prior.getExponent();
    }

    @Override
    public Parameter getSlabWidth() {
        return prior.getSlabWidth();
    }

    @Override
    public int getDimension() {
        return rates.getDimension();
    }
}

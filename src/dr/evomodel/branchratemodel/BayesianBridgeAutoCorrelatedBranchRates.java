/*
 * AutoCorrelatedBranchRatesDistribution.java
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

import dr.evolution.tree.MutableTreeModel;
import dr.evolution.tree.TreeTrait;
import dr.evomodel.tree.TreeParameterModel;
import dr.inference.distribution.shrinkage.BayesianBridgeDistributionModel;
import dr.inference.distribution.shrinkage.BayesianBridgeStatisticsProvider;
import dr.inference.model.*;

/**
 * @author Marc A. Suchard
 * @author Alexander Fisher
 */
public class BayesianBridgeAutoCorrelatedBranchRates extends AutoCorrelatedBranchRatesDistribution
        implements BayesianBridgeStatisticsProvider {

    private final BayesianBridgeDistributionModel distribution;
    private final TreeParameterModel treeParameterModel;

    public BayesianBridgeAutoCorrelatedBranchRates(String name,
                                                   DifferentiableBranchRates branchRateModel,
                                                   BayesianBridgeDistributionModel distribution,
                                                   BranchVarianceScaling scaling,
                                                   boolean takeLogBeforeIncrement, boolean operateOnIncrements) {
        super(name, branchRateModel, distribution, scaling, takeLogBeforeIncrement, operateOnIncrements);
        this.distribution = distribution;
        this.treeParameterModel = new TreeParameterModel(
                (MutableTreeModel)branchRateModel.getTree(), getLocalScale(),
                false, TreeTrait.Intent.NODE);
        addModel(treeParameterModel);
    }

    @Override
    public double getCoefficient(int i) {
        return getIncrement(i);
    }

    @Override
    public Parameter getGlobalScale() {
        return distribution.getGlobalScale();
    }

    @Override
    public Parameter getLocalScale() {
        return distribution.getLocalScale();
    }

    @Override
    public Parameter getExponent() {
        return distribution.getExponent();
    }

    @Override
    public Parameter getSlabWidth() {
        return distribution.getSlabWidth();
    }
}

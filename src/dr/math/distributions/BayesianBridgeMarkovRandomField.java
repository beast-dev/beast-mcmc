/*
 * BayesianBridgeMarkovRandomField.java
 *
 * Copyright (c) 2002-2023 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.math.distributions;

import dr.inference.distribution.RandomField;
import dr.inference.distribution.shrinkage.BayesianBridgeStatisticsProvider;
import dr.inference.distribution.shrinkage.JointBayesianBridgeDistributionModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;

/**
 * @author Marc Suchard
 * @author Yucai Shao
 * @author Andy Magee
 */
public class BayesianBridgeMarkovRandomField extends GaussianMarkovRandomField implements
        BayesianBridgeStatisticsProvider {

    public static final String TYPE = "BayesianBridgeMarkovRandomField";

    private final JointBayesianBridgeDistributionModel bayesianBridge;

    public BayesianBridgeMarkovRandomField(String name,
                                           JointBayesianBridgeDistributionModel bayesianBridge,
                                           Parameter mean,
                                           RandomField.WeightProvider weightProvider,
                                           boolean matchPseudoDeterminant) {
        super(name, bayesianBridge.getDimension(), null, mean, weightProvider, matchPseudoDeterminant);

        this.bayesianBridge = bayesianBridge;
        addModel(bayesianBridge);
    }

    @Override
    SymmetricTriDiagonalMatrix getQ() {
        if (!qKnown) {
//            bayesianBridge.getStandardDeviation(i);
            throw new RuntimeException("Not yet implemented");
        }
        return Q;
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == bayesianBridge) {
            qKnown = false;
        } else {
            throw new IllegalArgumentException("Unknown model");
        }
    }

    @Override
    public double getCoefficient(int i) {
        throw new RuntimeException("Should never call");
    }

    @Override
    public Parameter getGlobalScale() {
        return bayesianBridge.getGlobalScale();
    }

    @Override
    public Parameter getLocalScale() {
        return bayesianBridge.getLocalScale();
    }

    @Override
    public Parameter getExponent() {
        return bayesianBridge.getExponent();
    }

    @Override
    public Parameter getSlabWidth() {
        return bayesianBridge.getSlabWidth();
    }
}

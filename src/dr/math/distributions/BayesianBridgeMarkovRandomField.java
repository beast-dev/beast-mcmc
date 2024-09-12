/*
 * BayesianBridgeMarkovRandomField.java
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
        super(name, bayesianBridge.getDimension(), null, mean, null, weightProvider, matchPseudoDeterminant);

        this.bayesianBridge = bayesianBridge;
        addModel(bayesianBridge);
    }

    @Override
    protected SymmetricTriDiagonalMatrix getQ() {
        if (!qKnown) {

            final double[] diagonal = Q.diagonal;
            final double[] offDiagonal = Q.offDiagonal;

            double sd = bayesianBridge.getStandardDeviation(0);
            offDiagonal[0] = -1.0 / (sd * sd);
            diagonal[0] = -offDiagonal[0];

            for (int i = 1; i < dim - 1; ++i) {
                sd = bayesianBridge.getStandardDeviation(i);
                offDiagonal[i] = -1.0 / (sd * sd);
                diagonal[i] = -(offDiagonal[i - 1] + offDiagonal[i]);
            }

            diagonal[dim - 1] = -offDiagonal[dim - 2];

            // TODO Update for lambda != 1 and for weights

            qKnown = true;
        }
        return Q;
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == bayesianBridge) {
            qKnown = false;
            // TODO do we need a fireModelChangedEvent()?
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

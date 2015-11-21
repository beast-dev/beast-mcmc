/*
 * MultiEpochExponentialModel.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.coalescent;

import dr.evolution.coalescent.DemographicFunction;
import dr.evolution.coalescent.ExponentialExponential;
import dr.evolution.coalescent.MultiEpochExponential;
import dr.inference.model.Parameter;

/**
 * @author Marc A. Suchard
 * @author Andrew Rambaut
 */
public class MultiEpochExponentialModel extends DemographicModel {

    /**
     * Construct demographic model with default settings
     */
    public MultiEpochExponentialModel(String name,
                                      Parameter N0Parameter,
                                      Parameter growthRateParameter,
                                      Parameter transitionTimeParameter,
                                      Type units) {

        super(name);

        int numEpoch = growthRateParameter.getDimension();
        if (numEpoch != transitionTimeParameter.getDimension() + 1) {
            throw new IllegalArgumentException("Invalid parameter dimensions");
        }

        multiEpochExponential = new MultiEpochExponential(units, numEpoch);

        this.N0Parameter = N0Parameter;
        addVariable(N0Parameter);
        N0Parameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

        this.growthRateParameter = growthRateParameter;
        addVariable(growthRateParameter);
        growthRateParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY,
                growthRateParameter.getDimension()));

        this.transitionTimeParameter = transitionTimeParameter;
        addVariable(transitionTimeParameter);
        transitionTimeParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY,
                0.0, transitionTimeParameter.getDimension()));

        setUnits(units);
    }


    // general functions

    public DemographicFunction getDemographicFunction() {
        multiEpochExponential.setN0(N0Parameter.getParameterValue(0));

        for (int i = 0; i < growthRateParameter.getDimension(); ++i) {
            multiEpochExponential.setGrowthRate(i, growthRateParameter.getParameterValue(i));
        }

        double totalTime = 0.0;
        for (int i = 0; i < transitionTimeParameter.getDimension(); ++i) {
            totalTime += transitionTimeParameter.getParameterValue(i);
            multiEpochExponential.setTransitionTime(i, totalTime);
        }

        return multiEpochExponential;
    }

    //
    // protected stuff
    //

    Parameter N0Parameter = null;
    Parameter growthRateParameter = null;
    Parameter transitionTimeParameter = null;
    MultiEpochExponential multiEpochExponential = null;
}
/*
 * ExponentialExponentialModel.java
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
import dr.evomodelxml.coalescent.ExponentialExponentialModelParser;
import dr.inference.model.Parameter;

/**
 * Exponential growth followed by a different phase of exponential growth.
 *
 * @author Andrew Rambaut
 */
public class ExponentialExponentialModel extends DemographicModel {

    //
    // Public stuff
    //

    /**
     * Construct demographic model with default settings
     */
    public ExponentialExponentialModel(Parameter N0Parameter,
                                       Parameter growthRateParameter,
                                       Parameter ancestralGrowthRateParameter,
                                       Parameter transitionTimeParameter,
                                       Type units) {

        this(ExponentialExponentialModelParser.EXPONENTIAL_EXPONENTIAL_MODEL,
                N0Parameter,
                growthRateParameter,
                ancestralGrowthRateParameter,
                transitionTimeParameter,
                units);
    }

    /**
     * Construct demographic model with default settings
     */
    public ExponentialExponentialModel(String name, Parameter N0Parameter,
                                       Parameter growthRateParameter,
                                       Parameter ancestralGrowthRateParameter,
                                       Parameter transitionTimeParameter,
                                       Type units) {

        super(name);

        exponentialExponential = new ExponentialExponential(units);

        this.N0Parameter = N0Parameter;
        addVariable(N0Parameter);
        N0Parameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

        this.growthRateParameter = growthRateParameter;
        addVariable(growthRateParameter);
        growthRateParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, 1));

        this.ancestralGrowthRateParameter = ancestralGrowthRateParameter;
        addVariable(ancestralGrowthRateParameter);
        ancestralGrowthRateParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, 1));

        this.transitionTimeParameter = transitionTimeParameter;
        addVariable(transitionTimeParameter);
        transitionTimeParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY,
                Double.NEGATIVE_INFINITY, 1));

        setUnits(units);
    }


    // general functions

    public DemographicFunction getDemographicFunction() {
        exponentialExponential.setN0(N0Parameter.getParameterValue(0));

        exponentialExponential.setGrowthRate(growthRateParameter.getParameterValue(0));

        exponentialExponential.setAncestralGrowthRate(ancestralGrowthRateParameter.getParameterValue(0));

        exponentialExponential.setTransitionTime(transitionTimeParameter.getParameterValue(0));

        return exponentialExponential;
    }

    //
    // protected stuff
    //

    Parameter N0Parameter = null;
    Parameter growthRateParameter = null;
    Parameter ancestralGrowthRateParameter = null;
    Parameter transitionTimeParameter = null;
    ExponentialExponential exponentialExponential = null;
}
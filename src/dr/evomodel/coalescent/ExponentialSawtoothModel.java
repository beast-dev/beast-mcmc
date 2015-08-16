/*
 * ExponentialSawtoothModel.java
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
import dr.evolution.coalescent.ExponentialSawtooth;
import dr.evomodelxml.coalescent.ExponentialSawtoothModelParser;
import dr.inference.model.Parameter;

/**
 * @author Alexei Drummond
 * @version $Id: ExponentialSawtoothModel.java,v 1.4 2005/04/11 11:24:39 alexei Exp $
 */
public class ExponentialSawtoothModel extends DemographicModel {

    //
    // Public stuff
    //
    /**
     * Construct demographic model with default settings
     */
    public ExponentialSawtoothModel(Parameter N0Parameter, Parameter growthRateParameter, Parameter wavelengthParameter, Parameter offsetParameter, Type units) {

        this(ExponentialSawtoothModelParser.EXPONENTIAL_SAWTOOTH, N0Parameter, growthRateParameter, wavelengthParameter, offsetParameter, units);
    }

    /**
     * Construct demographic model with default settings
     */
    public ExponentialSawtoothModel(String name, Parameter N0Parameter, Parameter growthRateParameter, Parameter wavelengthParameter, Parameter offsetParameter, Type units) {

        super(name);

        expSaw = new ExponentialSawtooth(units);

        this.N0Parameter = N0Parameter;
        addVariable(N0Parameter);
        N0Parameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, Double.MIN_VALUE, 1));

        this.growthRateParameter = growthRateParameter;
        addVariable(growthRateParameter);
        growthRateParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, -Double.MAX_VALUE, 1));

        this.wavelengthParameter = wavelengthParameter;
        addVariable(wavelengthParameter);
        wavelengthParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, Double.MIN_VALUE, 1));

        this.offsetParameter = offsetParameter;
        addVariable(offsetParameter);
        offsetParameter.addBounds(new Parameter.DefaultBounds(1.0, -1.0, 1));

        setUnits(units);
    }


    // general functions

    public DemographicFunction getDemographicFunction() {
        expSaw.setN0(N0Parameter.getParameterValue(0));
        expSaw.setGrowthRate(growthRateParameter.getParameterValue(0));
        expSaw.setWavelength(wavelengthParameter.getParameterValue(0));

        double offset = offsetParameter.getParameterValue(0);
        if (offset < 0.0) {
            offset += 1.0;
        }
        expSaw.setOffset(offset);

        return expSaw;
    }
    //
    // protected stuff
    //

    Parameter N0Parameter = null;
    Parameter growthRateParameter = null;
    Parameter wavelengthParameter = null;
    Parameter offsetParameter = null;
    ExponentialSawtooth expSaw = null;
}

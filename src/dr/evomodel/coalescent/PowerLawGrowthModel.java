/*
 * PowerLawGrowthModel.java
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
import dr.evolution.coalescent.PowerLawGrowth;
import dr.evomodelxml.coalescent.PowerLawGrowthModelParser;
import dr.inference.model.Parameter;

/**
 * This class models an exponentially growing (or shrinking) population
 * (Parameters: N0=present-day population size; r=growth rate).
 * This model is nested with the constant-population size model (r=0).
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @version $Id: ExponentialGrowthModel.java,v 1.14 2005/05/24 20:25:57 rambaut Exp $
 */
public class PowerLawGrowthModel extends DemographicModel {

    //
    // Public stuff
    //
    /**
     * Construct demographic model with default settings
     */
    public PowerLawGrowthModel(Parameter N0Parameter, Parameter growthRateParameter,
                               Type units) {

        this(PowerLawGrowthModelParser.POWER_LAW_GROWTH_MODEL, N0Parameter, growthRateParameter, units);
    }

    /**
     * Construct demographic model with default settings
     */
    public PowerLawGrowthModel(String name, Parameter N0Parameter, Parameter powerParameter,
                               Type units) {

        super(name);

        powerLawGrowth = new PowerLawGrowth(units);

        this.N0Parameter = N0Parameter;
        addVariable(N0Parameter);
        N0Parameter.addBounds(new Parameter.DefaultBounds(Double.MAX_VALUE, 0.0, 1));

        this.powerParameter = powerParameter;
        addVariable(powerParameter);
        powerParameter.addBounds(new Parameter.DefaultBounds(Double.MAX_VALUE, 1, 1));

        setUnits(units);
    }


    // general functions

    public DemographicFunction getDemographicFunction() {
        powerLawGrowth.setN0(N0Parameter.getParameterValue(0));


        double r = powerParameter.getParameterValue(0);
        powerLawGrowth.setR(r);

        return powerLawGrowth;
    }
    //
    // protected stuff
    //

    Parameter N0Parameter = null;
    Parameter powerParameter = null;
    PowerLawGrowth powerLawGrowth = null;
}

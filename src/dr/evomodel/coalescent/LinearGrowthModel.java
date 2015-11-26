/*
 * LinearGrowthModel.java
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
import dr.evolution.coalescent.LinearGrowth;
import dr.evomodelxml.coalescent.LinearGrowthModelParser;
import dr.inference.model.Parameter;

/**
 * A wrapper for LinearGrowth.
 *
 * @author Matthew Hall
 *
 */
public class LinearGrowthModel extends DemographicModel {
    //
    // Public stuff
    //
    /**
     * Construct demographic model with default settings
     */
    public LinearGrowthModel(Parameter slopeParameter, Type units) {

        this(LinearGrowthModelParser.LINEAR_GROWTH_MODEL, slopeParameter, units);
    }

    /**
     * Construct demographic model with default settings
     */
    public LinearGrowthModel(String name, Parameter slopeParameter, Type units) {

        super(name);

        linearGrowth = new LinearGrowth(units);

        this.slopeParameter = slopeParameter;
        addVariable(slopeParameter);
        slopeParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
        setUnits(units);
    }

    // general functions

    public DemographicFunction getDemographicFunction() {
        linearGrowth.setN0(slopeParameter.getParameterValue(0));
        return linearGrowth;
    }

    //
    // protected stuff
    //

    private Parameter slopeParameter;
    private LinearGrowth linearGrowth = null;
}

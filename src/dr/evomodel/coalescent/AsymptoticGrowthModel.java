/*
 * AsymptoticGrowthModel.java
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
import dr.evolution.coalescent.FlexibleGrowth;
import dr.evomodelxml.coalescent.AsymptoticGrowthModelParser;
import dr.inference.model.Parameter;

/**
 * Growth starts at zero at time zero, peaks and declines
 *
 */
public class AsymptoticGrowthModel extends DemographicModel {

    //
    // Public stuff
    //

    /**
     * Construct demographic model with default settings
     */
    public AsymptoticGrowthModel(Parameter asymptoteValueParameter, Parameter shapeParameter, Type units) {

        this(AsymptoticGrowthModelParser.ASYMPTOTIC_GROWTH_MODEL, asymptoteValueParameter, shapeParameter, units);
    }

    /**
     * Construct demographic model with default settings
     */
    public AsymptoticGrowthModel(String name, Parameter asymptoteValueParameter, Parameter shapeParameter,
                                 Type units) {

        super(name);

        flexibleGrowth = new FlexibleGrowth(units);

        this.asyptoteValue = asymptoteValueParameter;
        addVariable(asymptoteValueParameter);
        asymptoteValueParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));


        this.shapeParameter = shapeParameter;
        addVariable(shapeParameter);
        shapeParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0, 1));


        setUnits(units);
    }


    // general functions

    public DemographicFunction getDemographicFunction() {

        double asymptoteValue = asyptoteValue.getParameterValue(0);
        double shapeValue = shapeParameter.getParameterValue(0);



        double flexibleN0 = asymptoteValue/shapeValue;


        flexibleGrowth.setN0(flexibleN0);
        flexibleGrowth.setK(shapeValue);
        flexibleGrowth.setR(0);



        return flexibleGrowth;
    }

    //
    // protected stuff
    //

    Parameter asyptoteValue = null;
    Parameter shapeParameter = null;
    FlexibleGrowth flexibleGrowth = null;

}

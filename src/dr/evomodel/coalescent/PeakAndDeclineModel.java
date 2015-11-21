/*
 * PeakAndDeclineModel.java
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
import dr.evomodelxml.coalescent.PeakAndDeclineModelParser;
import dr.inference.model.Parameter;

/**
 * Growth starts at zero at time zero, peaks and declines
 *
 */
public class PeakAndDeclineModel extends DemographicModel {

    //
    // Public stuff
    //

    /**
     * Construct demographic model with default settings
     */
    public PeakAndDeclineModel(Parameter peakValueParameter, Parameter shapeParameter, Parameter peakTimeParameter,
                               Type units) {

        this(PeakAndDeclineModelParser.PEAK_AND_DECLINE_MODEL, peakValueParameter, shapeParameter, peakTimeParameter,
                units);
    }

    /**
     * Construct demographic model with default settings
     */
    public PeakAndDeclineModel(String name, Parameter peakValueParameter, Parameter shapeParameter,
                               Parameter peakTimeParameter, Type units) {

        super(name);

        flexibleGrowth = new FlexibleGrowth(units);

        this.peakValueParameter = peakValueParameter;
        addVariable(peakValueParameter);
        peakValueParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

        this.peakTimeParameter = peakTimeParameter;
        addVariable(peakTimeParameter);
        peakTimeParameter.addBounds(new Parameter.DefaultBounds(0, Double.NEGATIVE_INFINITY, 1));

        this.shapeParameter = shapeParameter;
        addVariable(shapeParameter);
        shapeParameter.addBounds(new Parameter.DefaultBounds(0, Double.NEGATIVE_INFINITY, 1));


        setUnits(units);
    }


    // general functions

    public DemographicFunction getDemographicFunction() {

        double peakTimeValue = peakTimeParameter.getParameterValue(0);
        double peakValueValue = peakValueParameter.getParameterValue(0);
        double shapeValue = shapeParameter.getParameterValue(0);

        double flexibleN0 = peakValueValue*(1-shapeValue)/(shapeValue*peakTimeValue);
        double flexibleK = (-shapeValue/Math.pow(-peakTimeValue, shapeValue-1));


        flexibleGrowth.setN0(flexibleN0);
        flexibleGrowth.setK(flexibleK);
        flexibleGrowth.setR(shapeValue);



        return flexibleGrowth;
    }

    //
    // protected stuff
    //

    Parameter peakValueParameter = null;
    Parameter shapeParameter = null;
    Parameter peakTimeParameter = null;
    FlexibleGrowth flexibleGrowth = null;

}

/*
 * LogisticGrowthModel.java
 *
 * Copyright (C) 2002-2009 Alexei Drummond and Andrew Rambaut
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
 * BEAST is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.evomodel.coalescent;

import dr.evolution.coalescent.DemographicFunction;
import dr.evolution.coalescent.FlexibleGrowth;
import dr.evolution.coalescent.LogisticGrowth;
import dr.evomodelxml.coalescent.LogisticGrowthModelParser;
import dr.inference.model.Parameter;

/**
 * Logistic growth.
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @version $Id: LogisticGrowthModel.java,v 1.21 2005/05/24 20:25:57 rambaut Exp $
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

        this(LogisticGrowthModelParser.LOGISTIC_GROWTH_MODEL, peakValueParameter, shapeParameter, peakTimeParameter, units);
    }

    /**
     * Construct demographic model with default settings
     */
    public PeakAndDeclineModel(String name, Parameter peakValueParameter, Parameter shapeParameter, Parameter peakTimeParameter,
                               Type units) {

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

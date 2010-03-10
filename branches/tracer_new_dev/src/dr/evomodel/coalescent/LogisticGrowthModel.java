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
public class LogisticGrowthModel extends DemographicModel {

    //
    // Public stuff
    //

    /**
     * Construct demographic model with default settings
     */
    public LogisticGrowthModel(Parameter N0Parameter, Parameter growthRateParameter,
                               Parameter shapeParameter, double alpha, Type units,
                               boolean usingGrowthRate) {

        this(LogisticGrowthModelParser.LOGISTIC_GROWTH_MODEL, N0Parameter, growthRateParameter, shapeParameter, alpha, units, usingGrowthRate);
    }

    /**
     * Construct demographic model with default settings
     */
    public LogisticGrowthModel(String name, Parameter N0Parameter, Parameter growthRateParameter, Parameter shapeParameter, double alpha, Type units, boolean usingGrowthRate) {

        super(name);

        logisticGrowth = new LogisticGrowth(units);

        this.N0Parameter = N0Parameter;
        addVariable(N0Parameter);
        N0Parameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

        this.growthRateParameter = growthRateParameter;
        addVariable(growthRateParameter);
        growthRateParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

        this.shapeParameter = shapeParameter;
        addVariable(shapeParameter);
        shapeParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

        this.alpha = alpha;
        this.usingGrowthRate = usingGrowthRate;

        setUnits(units);
    }


    // general functions

    public DemographicFunction getDemographicFunction() {

        logisticGrowth.setN0(N0Parameter.getParameterValue(0));

        if (usingGrowthRate) {
            double r = growthRateParameter.getParameterValue(0);
            logisticGrowth.setGrowthRate(r);
        } else {
            double doublingTime = growthRateParameter.getParameterValue(0);
            logisticGrowth.setDoublingTime(doublingTime);
        }

        logisticGrowth.setTime50(shapeParameter.getParameterValue(0));

        return logisticGrowth;
    }

    //
    // protected stuff
    //

    Parameter N0Parameter = null;
    Parameter growthRateParameter = null;
    Parameter shapeParameter = null;
    double alpha = 0.5;
    LogisticGrowth logisticGrowth = null;
    boolean usingGrowthRate = true;
}

/*
 * ExponentialGrowthModelParser.java
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

package dr.evomodelxml.coalescent;

import dr.evolution.util.Units;
import dr.evomodel.coalescent.ExponentialPopulationSizeModel;
import dr.evomodel.coalescent.demographicmodel.ExponentialGrowthModel;
import dr.evoxml.util.XMLUnits;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * Parses an element from an DOM document into a ExponentialGrowth.
 */
public class ExponentialPopulationSizeModelParser extends AbstractXMLObjectParser {

    public static String EXPONENTIAL_POPULATION_SIZE_MODEL = "exponentialPopulationSize";
    public static String LOG_POPULATION_SIZE = "logPopulationSize";

    public static String GROWTH_RATE = "growthRate";
    public static String DOUBLING_TIME = "doublingTime";


    public String getParserName() {
        return EXPONENTIAL_POPULATION_SIZE_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Units.Type units = XMLUnits.Utils.getUnitsAttr(xo);

        Parameter N0Param = null;

        if (xo.hasChildNamed(LOG_POPULATION_SIZE)) {
            N0Param = (Parameter)xo.getElementFirstChild(LOG_POPULATION_SIZE);
        }

        Parameter rateParameter = (Parameter) xo.getElementFirstChild(GROWTH_RATE);

        // ExponentialPopulationSizeModel provides a model with N0 in log space
        return new ExponentialPopulationSizeModel(N0Param, rateParameter, units);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A demographic model of exponential growth.";
    }

    public Class getReturnType() {
        return ExponentialGrowthModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(LOG_POPULATION_SIZE,  Parameter.class, "the log population size parameter", true),
            new ElementRule(GROWTH_RATE,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)},
                    "A value of zero represents a constant population size, negative values represent decline towards the present, " +
                            "positive numbers represents exponential growth towards the present. " +
                            "A random walk operator is recommended for this parameter with a starting value of 0.0 and no upper or lower limits."),
            XMLUnits.SYNTAX_RULES[0]
    };


}

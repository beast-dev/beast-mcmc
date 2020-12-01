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

package dr.evomodelxml.coalescent.demographicmodel;

import dr.evolution.util.Units;
import dr.evomodel.coalescent.ExponentialPopulationSizeModel;
import dr.evomodel.coalescent.demographicmodel.ExponentialGrowthModel;
import dr.evoxml.util.XMLUnits;
import dr.inference.model.Parameter;
import dr.xml.*;

import static dr.evomodelxml.coalescent.demographicmodel.ConstantPopulationModelParser.LOG_SPACE;

/**
 * Parses an element from an DOM document into a ExponentialGrowth.
 */
public class ExponentialGrowthModelParser extends AbstractXMLObjectParser {

    public static String POPULATION_SIZE = "populationSize";
    public static String EXPONENTIAL_GROWTH_MODEL = "exponentialGrowth";

    public static String GROWTH_RATE = "growthRate";
    public static String DOUBLING_TIME = "doublingTime";


    public String getParserName() {
        return EXPONENTIAL_GROWTH_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Units.Type units = XMLUnits.Utils.getUnitsAttr(xo);

        XMLObject cxo = xo.getChild(POPULATION_SIZE);
        boolean logSpace = cxo.getAttribute(LOG_SPACE, false);
        Parameter N0Param = (Parameter) cxo.getChild(Parameter.class);

        Parameter rParam;
        boolean usingGrowthRate = true;

        if (xo.getChild(GROWTH_RATE) != null) {
            cxo = xo.getChild(GROWTH_RATE);
            rParam = (Parameter) cxo.getChild(Parameter.class);
        } else {
            cxo = xo.getChild(DOUBLING_TIME);
            rParam = (Parameter) cxo.getChild(Parameter.class);
            usingGrowthRate = false;
        }

        if (logSpace) {
            if (!usingGrowthRate) {
                throw new XMLParseException("Doubling time parameterization is not compatible with log population size");
            }
            // ExponentialPopulationSizeModel provides a model with N0 in log space
            return new ExponentialPopulationSizeModel(N0Param, rParam, units);
        } else {
            return new ExponentialGrowthModel(N0Param, rParam, units, usingGrowthRate);
        }
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
            new ElementRule(POPULATION_SIZE,
                    new XMLSyntaxRule[] {
                            AttributeRule.newBooleanRule(LOG_SPACE, true, "Is this parameter in log space?"),
                            new ElementRule(Parameter.class)}),
            new XORRule(

                    new ElementRule(GROWTH_RATE,
                            new XMLSyntaxRule[]{new ElementRule(Parameter.class)},
                            "A value of zero represents a constant population size, negative values represent decline towards the present, " +
                                    "positive numbers represents exponential growth towards the present. " +
                                    "A random walk operator is recommended for this parameter with a starting value of 0.0 and no upper or lower limits."),
                    new ElementRule(DOUBLING_TIME,
                            new XMLSyntaxRule[]{new ElementRule(Parameter.class)},
                            "This parameter determines the doubling time.")
            ),
            XMLUnits.SYNTAX_RULES[0]
    };


}

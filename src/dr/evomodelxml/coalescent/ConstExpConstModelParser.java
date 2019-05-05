/*
 * ConstExpConstModelParser.java
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
import dr.evomodel.coalescent.demographicmodels.ConstExpConstModel;
import dr.evomodel.coalescent.demographicmodels.ConstantExponentialModel;
import dr.evoxml.util.XMLUnits;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * Parses an element from an DOM document into a ConstantExponentialModel.
 *
 * There are 2 parameterization depending on whether a growth rate parameter
 * or an ancestral population size parameter are specified.
 */
public class ConstExpConstModelParser extends AbstractXMLObjectParser {

    public static String CONST_EXP_CONST_MODEL = "constExpConst";
    public static String POPULATION_SIZE = "populationSize";
    public static String GROWTH_RATE = "growthRate";
    public static String ANCESTRAL_POPULATION_SIZE = "ancestralPopulationSize";
    public static String FINAL_PHASE_START_TIME = "finalPhaseStartTime";
    public static String GROWTH_PHASE_TIME = "growthPhaseTime";
    public static String USE_NUMERICAL_INTEGRATION = "useNumericalIntegration";

    public String getParserName() {
        return CONST_EXP_CONST_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Units.Type units = XMLUnits.Utils.getUnitsAttr(xo);

        XMLObject cxo = xo.getChild(POPULATION_SIZE);
        Parameter N0Param = (Parameter) cxo.getChild(Parameter.class);

        Parameter N1Param = null;
        Parameter growthRateParam = null;

        if (xo.hasChildNamed(ANCESTRAL_POPULATION_SIZE)) {
            cxo = xo.getChild(ANCESTRAL_POPULATION_SIZE);
            N1Param = (Parameter) cxo.getChild(Parameter.class);
        } else {
            cxo = xo.getChild(GROWTH_RATE);
            growthRateParam = (Parameter) cxo.getChild(Parameter.class);
        }

        cxo = xo.getChild(FINAL_PHASE_START_TIME);
        Parameter timeParam = (Parameter) cxo.getChild(Parameter.class);

        cxo = xo.getChild(GROWTH_PHASE_TIME);
        Parameter epochParam = (Parameter) cxo.getChild(Parameter.class);

        boolean useNumericalIntegrator = false;
        if (xo.hasAttribute(USE_NUMERICAL_INTEGRATION)) {
            useNumericalIntegrator = xo.getBooleanAttribute(USE_NUMERICAL_INTEGRATION);
        }

        return new ConstExpConstModel(N0Param, N1Param, growthRateParam, timeParam, epochParam, useNumericalIntegrator, units);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A demographic model of constant population size followed by exponential growth.";
    }

    public Class getReturnType() {
        return ConstantExponentialModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            XMLUnits.SYNTAX_RULES[0],
            AttributeRule.newBooleanRule(USE_NUMERICAL_INTEGRATION, true),
            new ElementRule(POPULATION_SIZE,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new XORRule(
                    new ElementRule(ANCESTRAL_POPULATION_SIZE,
                            new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
                    new ElementRule(GROWTH_RATE,
                            new XMLSyntaxRule[]{new ElementRule(Parameter.class)})
            ),
            new ElementRule(FINAL_PHASE_START_TIME,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(GROWTH_PHASE_TIME,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
    };


}

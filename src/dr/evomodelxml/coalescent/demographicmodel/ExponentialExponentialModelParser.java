/*
 * ExponentialExponentialModelParser.java
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
import dr.evomodel.coalescent.demographicmodel.ExponentialExponentialModel;
import dr.evoxml.util.XMLUnits;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * Parses an element from an DOM document into a ExponentialExponentialModel.
 */
public class ExponentialExponentialModelParser extends AbstractXMLObjectParser {

    public static final String EXPONENTIAL_EXPONENTIAL_MODEL = "exponentialExponential";
    public static final String POPULATION_SIZE = "populationSize";
    public static final String TRANSITION_POPULATION_SIZE = "transitionPopulationSize";
    public static final String TRANSITION_TIME = "transitionTime";
    public static final String ANCESTRAL_GROWTH_RATE = "ancestralGrowthRate";
    public static final String GROWTH_RATE = "growthRate";

    public String getParserName() {
        return EXPONENTIAL_EXPONENTIAL_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Units.Type units = XMLUnits.Utils.getUnitsAttr(xo);

        Parameter N0Param = null;
        Parameter N1Param = null;

        if (xo.hasChildNamed(POPULATION_SIZE)) {
            XMLObject cxo = xo.getChild(POPULATION_SIZE);
            N0Param = (Parameter) cxo.getChild(Parameter.class);
        } else {
            // must be the TRANSITION_POPULATION_SIZE
            XMLObject cxo = xo.getChild(TRANSITION_POPULATION_SIZE);
            N1Param = (Parameter) cxo.getChild(Parameter.class);
        }

        XMLObject cxo = xo.getChild(GROWTH_RATE);
        Parameter growthParam = (Parameter) cxo.getChild(Parameter.class);

        cxo = xo.getChild(ANCESTRAL_GROWTH_RATE);
        Parameter ancestralGrowthParam = (Parameter) cxo.getChild(Parameter.class);

        cxo = xo.getChild(TRANSITION_TIME);
        Parameter timeParam = (Parameter) cxo.getChild(Parameter.class);

        return new ExponentialExponentialModel(N0Param, N1Param, growthParam, ancestralGrowthParam, timeParam, units);
    }

//************************************************************************
// AbstractXMLObjectParser implementation
//************************************************************************

    public String getParserDescription() {
        return "A demographic model of exponential growth followed by a different exponential growth.";
    }

    public Class getReturnType() {
        return ExponentialExponentialModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            XMLUnits.SYNTAX_RULES[0],
            new XORRule(
                    new ElementRule(POPULATION_SIZE,
                            new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
                    new ElementRule(TRANSITION_POPULATION_SIZE,
                            new XMLSyntaxRule[]{new ElementRule(Parameter.class)})
            ),
            new ElementRule(GROWTH_RATE,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(ANCESTRAL_GROWTH_RATE,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(TRANSITION_TIME,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
    };
}

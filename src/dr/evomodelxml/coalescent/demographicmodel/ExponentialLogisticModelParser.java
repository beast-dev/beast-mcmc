/*
 * ExponentialLogisticModelParser.java
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
import dr.evomodel.coalescent.demographicmodel.ExponentialLogisticModel;
import dr.evoxml.util.XMLUnits;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * Parses an element from an DOM document into a ExponentialGrowth.
 */
public class ExponentialLogisticModelParser extends AbstractXMLObjectParser {

    public static final String EXPONENTIAL_LOGISTIC_MODEL = "exponentialLogistic";
    public static final String POPULATION_SIZE = "populationSize";
    public static final String TRANSITION_TIME = "transitionTime";

    public static final String LOGISTIC_GROWTH_RATE = "logisticGrowthRate";
    public static final String LOGISTIC_SHAPE = "logisticShape";
    public static final String EXPONENTIAL_GROWTH_RATE = "exponentialGrowthRate";

    public static final String ALPHA = "alpha";

    public String getParserName() {
        return EXPONENTIAL_LOGISTIC_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Units.Type units = XMLUnits.Utils.getUnitsAttr(xo);

        XMLObject cxo = xo.getChild(POPULATION_SIZE);
        Parameter N0Param = (Parameter) cxo.getChild(Parameter.class);

        cxo = xo.getChild(LOGISTIC_GROWTH_RATE);
        Parameter logGrowthParam = (Parameter) cxo.getChild(Parameter.class);

        cxo = xo.getChild(LOGISTIC_SHAPE);
        Parameter shapeParam = (Parameter) cxo.getChild(Parameter.class);

        cxo = xo.getChild(EXPONENTIAL_GROWTH_RATE);
        Parameter expGrowthParam = (Parameter) cxo.getChild(Parameter.class);

        cxo = xo.getChild(TRANSITION_TIME);
        Parameter timeParam = (Parameter) cxo.getChild(Parameter.class);

        double alpha = xo.getDoubleAttribute(ALPHA);

        return new ExponentialLogisticModel(N0Param, logGrowthParam, shapeParam, expGrowthParam, timeParam, alpha, units);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A demographic model of exponential growth followed by logistic growth.";
    }

    public Class getReturnType() {
        return ExponentialLogisticModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            XMLUnits.SYNTAX_RULES[0],
            new ElementRule(POPULATION_SIZE,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(LOGISTIC_GROWTH_RATE,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(LOGISTIC_SHAPE,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(EXPONENTIAL_GROWTH_RATE,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(TRANSITION_TIME,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            AttributeRule.newDoubleRule(ALPHA)
    };
}

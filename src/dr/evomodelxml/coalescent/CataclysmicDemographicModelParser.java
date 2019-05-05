/*
 * CataclysmicDemographicModelParser.java
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
import dr.evomodel.coalescent.demographicmodels.CataclysmicDemographicModel;
import dr.evoxml.util.XMLUnits;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * Parses an element from an DOM document into a ExponentialGrowth.
 */
public class CataclysmicDemographicModelParser extends AbstractXMLObjectParser {

    public static final String POPULATION_SIZE = "populationSize";
    public static final String GROWTH_RATE = "growthRate";
    public static final String SPIKE_SIZE = "spikeFactor";
    public static final String TIME_OF_CATACLYSM = "timeOfCataclysm";
    public static final String DECLINE_RATE = "declineRate";

    public static final String CATACLYSM_MODEL = "cataclysm";

    public String getParserName() {
        return CATACLYSM_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Units.Type units = XMLUnits.Utils.getUnitsAttr(xo);

        XMLObject cxo = xo.getChild(POPULATION_SIZE);
        Parameter N0Param = (Parameter) cxo.getChild(Parameter.class);

        cxo = xo.getChild(GROWTH_RATE);
        Parameter rParam = (Parameter) cxo.getChild(Parameter.class);

        Parameter secondParam = null;
        boolean useSpike = true;

        if (xo.hasChildNamed(SPIKE_SIZE)) {
            cxo = xo.getChild(SPIKE_SIZE);
            secondParam = (Parameter) cxo.getChild(Parameter.class);
        } else if (xo.hasChildNamed(DECLINE_RATE)) {
            cxo = xo.getChild(DECLINE_RATE);
            secondParam = (Parameter) cxo.getChild(Parameter.class);
            useSpike = false;
        } else {
            throw new XMLParseException("Must provide either a spike factor or decline rate");
        }

        cxo = xo.getChild(TIME_OF_CATACLYSM);
        Parameter tParam = (Parameter) cxo.getChild(Parameter.class);

        return new CataclysmicDemographicModel(N0Param, secondParam, rParam, tParam, units, useSpike);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A demographic model of exponential growth.";
    }

    public Class getReturnType() {
        return CataclysmicDemographicModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(POPULATION_SIZE,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(GROWTH_RATE,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)},
                    "The rate of exponential growth before the cataclysmic event."),
            new XORRule(
                    new ElementRule(SPIKE_SIZE,
                        new XMLSyntaxRule[]{new ElementRule(Parameter.class)},
                        "The factor larger the population size was at its height."),
                    new ElementRule(DECLINE_RATE,
                            new XMLSyntaxRule[] { new ElementRule(Parameter.class)})),
            new ElementRule(TIME_OF_CATACLYSM,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)},
                    "The time of the cataclysmic event that lead to exponential decline."),
            XMLUnits.SYNTAX_RULES[0]
    };
}

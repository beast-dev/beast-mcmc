/*
 * ExpConstExpDemographicModelParser.java
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
import dr.evomodel.coalescent.demographicmodel.ExpConstExpDemographicModel;
import dr.evoxml.util.XMLUnits;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * Parses an element from an DOM document into a ExponentialGrowth.
 */
public class ExpConstExpDemographicModelParser extends AbstractXMLObjectParser {

    public static final String POPULATION_SIZE = "populationSize";
    public static final String RELATIVE_PLATEAU_SIZE = "relativePlateauSize";
    public static final String RELATIVE_TIME_OF_MODERN_GROWTH = "relTimeModGrowth";
    public static final String TIME_PLATEAU = "plateauStartTime";
    public static final String ANCIENT_GROWTH_RATE = "ancientGrowthRate";

    public static final String EXP_CONST_EXP_MODEL = "expConstExp";

    public String getParserName() {
        return EXP_CONST_EXP_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Units.Type units = XMLUnits.Utils.getUnitsAttr(xo);

        XMLObject cxo = xo.getChild(POPULATION_SIZE);
        Parameter N0Param = (Parameter) cxo.getChild(Parameter.class);

        cxo = xo.getChild(RELATIVE_PLATEAU_SIZE);
        Parameter N1Param = (Parameter) cxo.getChild(Parameter.class);

        cxo = xo.getChild(RELATIVE_TIME_OF_MODERN_GROWTH);
        Parameter rtParam = (Parameter) cxo.getChild(Parameter.class);

        cxo = xo.getChild(TIME_PLATEAU);
        Parameter tParam = (Parameter) cxo.getChild(Parameter.class);

        cxo = xo.getChild(ANCIENT_GROWTH_RATE);
        Parameter rParam = (Parameter) cxo.getChild(Parameter.class);

        return new ExpConstExpDemographicModel(N0Param, N1Param, rParam, tParam, rtParam, units);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A demographic model of exponential growth.";
    }

    public Class getReturnType() {
        return ExpConstExpDemographicModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(POPULATION_SIZE,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(RELATIVE_PLATEAU_SIZE,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)},
                    "The size of plateau relative to modern population size."),
            new ElementRule(RELATIVE_TIME_OF_MODERN_GROWTH,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)},
                    "The time spanned by modern growth phase relative to time back to start of plateau phase."),
            new ElementRule(TIME_PLATEAU,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)},
                    "The time of the start of plateauPhase."),
            new ElementRule(ANCIENT_GROWTH_RATE,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)},
                    "The growth rate of early growth phase"),
            XMLUnits.SYNTAX_RULES[0]
    };

}

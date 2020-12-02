/*
 * SIRModelParser.java
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

package dr.evomodel.epidemiology;

import dr.evolution.util.Units;
import dr.evoxml.util.XMLUnits;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * Parses an element from an DOM document into a SIR model.
 *
 * @author Trevor Bedford
 * @author Andrew Rambaut
 */
public class SIRModelParser extends AbstractXMLObjectParser {

    public static String SIR_MODEL = "sirEpidemiology";
    public static String REPRODUCTIVE_NUMBER = "reproductiveNumber";
    public static String RECOVERY_RATE = "recoveryRate";
    public static String HOST_POPULATION_SIZE = "hostPopulationSize";
    public static String PROPORTIONS = "proportions";

    public String getParserName() {
        return SIR_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Units.Type units = XMLUnits.Utils.getUnitsAttr(xo);

        XMLObject cxo = xo.getChild(REPRODUCTIVE_NUMBER);
        Parameter reproductiveNumberParameter = (Parameter) cxo.getChild(Parameter.class);

        cxo = xo.getChild(RECOVERY_RATE);
        Parameter recoveryRateParameter = (Parameter) cxo.getChild(Parameter.class);

        cxo = xo.getChild(HOST_POPULATION_SIZE);
        Parameter hostPopulationSizeParameter = (Parameter) cxo.getChild(Parameter.class);

        cxo = xo.getChild(PROPORTIONS);
        Parameter proportionsParameter = (Parameter) cxo.getChild(Parameter.class);

        return new SIRModel(reproductiveNumberParameter, recoveryRateParameter,
                hostPopulationSizeParameter, proportionsParameter, units);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A demographic model of epidemic spread.";
    }

    public Class getReturnType() {
        return SIRModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(REPRODUCTIVE_NUMBER,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(RECOVERY_RATE,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(HOST_POPULATION_SIZE,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(PROPORTIONS,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            XMLUnits.SYNTAX_RULES[0]
    };


}
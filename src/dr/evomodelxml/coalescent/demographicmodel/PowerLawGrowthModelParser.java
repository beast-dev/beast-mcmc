/*
 * PowerLawGrowthModelParser.java
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
import dr.evomodel.coalescent.demographicmodel.PowerLawGrowthModel;
import dr.evoxml.util.XMLUnits;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * Parses an element from an DOM document into a ExponentialGrowth.
 */
public class PowerLawGrowthModelParser extends AbstractXMLObjectParser {

    public static String N0 = "n0";
    public static String POWER_LAW_GROWTH_MODEL = "powerLawGrowth";

    public static String POWER = "power";


    public String getParserName() {
        return POWER_LAW_GROWTH_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Units.Type units = XMLUnits.Utils.getUnitsAttr(xo);

        XMLObject cxo = xo.getChild(N0);
        Parameter N0Param = (Parameter) cxo.getChild(Parameter.class);
        Parameter rParam;


        cxo = xo.getChild(POWER);
        rParam = (Parameter) cxo.getChild(Parameter.class);
        return new PowerLawGrowthModel(N0Param, rParam, units);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A demographic model of growth according to a power law.";
    }

    public Class getReturnType() {
        return PowerLawGrowthModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(N0,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),


            new ElementRule(POWER,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),

            XMLUnits.SYNTAX_RULES[0]
    };


}

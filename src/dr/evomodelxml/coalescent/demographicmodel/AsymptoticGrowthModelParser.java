/*
 * AsymptoticGrowthModelParser.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.evomodelxml.coalescent.demographicmodel;

import dr.evolution.util.Units;
import dr.evomodel.coalescent.demographicmodel.AsymptoticGrowthModel;
import dr.evomodel.coalescent.demographicmodel.LogisticGrowthModel;
import dr.evoxml.util.XMLUnits;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * Parses an element from an XMLObject into LogisticGrowthModel.
 */
public class AsymptoticGrowthModelParser extends AbstractXMLObjectParser {

    public static String ASYMPTOTE_VALUE = "asymptoteValue";
    public static String ASYMPTOTIC_GROWTH_MODEL = "asymptoticGrowth";

    public static String SHAPE = "shape";

    public String getParserName() {
        return ASYMPTOTIC_GROWTH_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Units.Type units = XMLUnits.Utils.getUnitsAttr(xo);

        XMLObject cxo = xo.getChild(ASYMPTOTE_VALUE);
        Parameter asymptoticValueParam = (Parameter) cxo.getChild(Parameter.class);


        Parameter rParam;

        cxo = xo.getChild(SHAPE);
        Parameter shapeParam = (Parameter) cxo.getChild(Parameter.class);


        return new AsymptoticGrowthModel(asymptoticValueParam, shapeParam, units);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "Logistic growth demographic model.";
    }

    public Class getReturnType() {
        return LogisticGrowthModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            XMLUnits.SYNTAX_RULES[0],
            new ElementRule(ASYMPTOTE_VALUE,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),


                    new ElementRule(SHAPE,
                            new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
    };
}

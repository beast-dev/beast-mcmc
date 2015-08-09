/*
 * TaxonRichnessBirthDeathModelParser.java
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

package dr.evomodelxml.speciation;

import dr.evolution.util.Units;
import dr.evomodel.speciation.TaxonRichnessBirthDeathModel;
import dr.evoxml.util.XMLUnits;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * @author Marc A. Suchard
 */
public class TaxonRichnessBirthDeathModelParser extends AbstractXMLObjectParser {

    public static final String TAXON_RICHNESS_MODEL = "randomLocalYuleModel";
    public static final String MEAN_RATE = "meanRate";
    public static final String BIRTH_RATE = "birthRates";
    public static final String BIRTH_RATE_INDICATORS = "indicators";
    public static final String RATES_AS_MULTIPLIERS = "ratesAsMultipliers";

    public String getParserName() {
        return TAXON_RICHNESS_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Units.Type units = XMLUnits.Utils.getUnitsAttr(xo);

        XMLObject cxo = xo.getChild(BIRTH_RATE);
        Parameter brParameter = (Parameter) cxo.getChild(Parameter.class);

        cxo = xo.getChild(BIRTH_RATE_INDICATORS);
        Parameter indicatorsParameter = (Parameter) cxo.getChild(Parameter.class);

        Parameter meanRate = (Parameter) xo.getElementFirstChild(MEAN_RATE);

        boolean ratesAsMultipliers = xo.getBooleanAttribute(RATES_AS_MULTIPLIERS);

        int dp = xo.getAttribute("dp", 4);

        return new TaxonRichnessBirthDeathModel(brParameter, indicatorsParameter, meanRate, ratesAsMultipliers, units, dp);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A speciation model of a Yule process whose rate can change at random nodes in the tree.";
    }

    public Class getReturnType() {
        return TaxonRichnessBirthDeathModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(BIRTH_RATE,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(BIRTH_RATE_INDICATORS,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(MEAN_RATE,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            AttributeRule.newBooleanRule(RATES_AS_MULTIPLIERS),
            AttributeRule.newIntegerRule("dp", true),
            XMLUnits.SYNTAX_RULES[0]
    };
}

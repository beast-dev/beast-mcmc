/*
 * GTRParser.java
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

package dr.oldevomodelxml.substmodel;

import dr.oldevomodel.substmodel.FrequencyModel;
import dr.oldevomodel.substmodel.GTR;
import dr.oldevomodel.substmodel.SubstitutionModel;
import dr.inference.model.Variable;
import dr.xml.*;

/**
 * Parses an element from an DOM document into a DemographicModel. Recognises
 * ConstantPopulation and ExponentialGrowth.
 */
public class GTRParser extends AbstractXMLObjectParser {
    public static final String GTR_MODEL = "gtrModel";

    public static final String A_TO_C = "rateAC";
    public static final String A_TO_G = "rateAG";
    public static final String A_TO_T = "rateAT";
    public static final String C_TO_G = "rateCG";
    public static final String C_TO_T = "rateCT";
    public static final String G_TO_T = "rateGT";

    public static final String FREQUENCIES = "frequencies";

    public String getParserName() {
        return GTR_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        XMLObject cxo = xo.getChild(FREQUENCIES);
        FrequencyModel freqModel = (FrequencyModel) cxo.getChild(FrequencyModel.class);

        Variable rateACValue = null;
        if (xo.hasChildNamed(A_TO_C)) {
            rateACValue = (Variable) xo.getElementFirstChild(A_TO_C);
        }
        Variable rateAGValue = null;
        if (xo.hasChildNamed(A_TO_G)) {
            rateAGValue = (Variable) xo.getElementFirstChild(A_TO_G);
        }
        Variable rateATValue = null;
        if (xo.hasChildNamed(A_TO_T)) {
            rateATValue = (Variable) xo.getElementFirstChild(A_TO_T);
        }
        Variable rateCGValue = null;
        if (xo.hasChildNamed(C_TO_G)) {
            rateCGValue = (Variable) xo.getElementFirstChild(C_TO_G);
        }
        Variable rateCTValue = null;
        if (xo.hasChildNamed(C_TO_T)) {
            rateCTValue = (Variable) xo.getElementFirstChild(C_TO_T);
        }
        Variable rateGTValue = null;
        if (xo.hasChildNamed(G_TO_T)) {
            rateGTValue = (Variable) xo.getElementFirstChild(G_TO_T);
        }
        int countNull = 0;
        if (rateACValue == null) countNull++;
        if (rateAGValue == null) countNull++;
        if (rateATValue == null) countNull++;
        if (rateCGValue == null) countNull++;
        if (rateCTValue == null) countNull++;
        if (rateGTValue == null) countNull++;

        if (countNull != 1)
            throw new XMLParseException("Only five parameters may be specified in GTR, leave exactly one out, the others will be specifed relative to the one left out.");
        return new GTR(rateACValue, rateAGValue, rateATValue, rateCGValue, rateCTValue, rateGTValue, freqModel);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A general reversible model of nucleotide sequence substitution.";
    }

    public String getExample() {

        return
                "<!-- A general time reversible model for DNA.                                          -->\n" +
                        "<!-- This element must have parameters for exactly five of the six rates               -->\n" +
                        "<!-- The sixth rate has an implied value of 1.0 and all other rates are relative to it -->\n" +
                        "<!-- This example parameterizes the rate matrix relative to the A<->G transition       -->\n" +
                        "<" + getParserName() + " id=\"gtr1\">\n" +
                        "	<" + FREQUENCIES + "> <frequencyModel idref=\"freqs\"/> </" + FREQUENCIES + ">\n" +
                        "	<" + A_TO_C + "> <parameter id=\"rateAC\" value=\"1.0\"/> </" + A_TO_C + ">\n" +
                        "	<" + A_TO_T + "> <parameter id=\"rateAT\" value=\"1.0\"/> </" + A_TO_T + ">\n" +
                        "	<" + C_TO_G + "> <parameter id=\"rateCG\" value=\"1.0\"/> </" + C_TO_G + ">\n" +
                        "	<" + C_TO_T + "> <parameter id=\"rateCT\" value=\"1.0\"/> </" + C_TO_T + ">\n" +
                        "	<" + G_TO_T + "> <parameter id=\"rateGT\" value=\"1.0\"/> </" + G_TO_T + ">\n" +
                        "</" + getParserName() + ">\n";
    }

    public Class getReturnType() {
        return SubstitutionModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(FREQUENCIES,
                    new XMLSyntaxRule[]{new ElementRule(FrequencyModel.class)}),
            new ElementRule(A_TO_C,
                    new XMLSyntaxRule[]{new ElementRule(Variable.class)}, true),
            new ElementRule(A_TO_G,
                    new XMLSyntaxRule[]{new ElementRule(Variable.class)}, true),
            new ElementRule(A_TO_T,
                    new XMLSyntaxRule[]{new ElementRule(Variable.class)}, true),
            new ElementRule(C_TO_G,
                    new XMLSyntaxRule[]{new ElementRule(Variable.class)}, true),
            new ElementRule(C_TO_T,
                    new XMLSyntaxRule[]{new ElementRule(Variable.class)}, true),
            new ElementRule(G_TO_T,
                    new XMLSyntaxRule[]{new ElementRule(Variable.class)}, true)
    };


}

/*
 * GTRParser.java
 *
 * Copyright (c) 2002-2016 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodelxml.substmodel;

import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.nucleotide.GTR;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.xml.*;

import static dr.evomodel.substmodel.nucleotide.GTR.*;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 */
public class GTRParser extends AbstractXMLObjectParser {
    public static final String GTR_MODEL = "gtrModel";

    public static final String FREQUENCIES = "frequencies";
    public static final String RATES = "rates";

    public static final String A_TO_C = "rateAC";
    public static final String A_TO_G = "rateAG";
    public static final String A_TO_T = "rateAT";
    public static final String C_TO_G = "rateCG";
    public static final String C_TO_T = "rateCT";
    public static final String G_TO_T = "rateGT";

    public String getParserName() {
        return GTR_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        FrequencyModel freqModel = (FrequencyModel) xo.getElementFirstChild(FREQUENCIES);

        Parameter rates = null;

        if (xo.hasChildNamed(RATES)) {
            rates = (Parameter) xo.getElementFirstChild(RATES);
            rates.setDimensionNames(new String[] {
                    rates.getId() + "." + A_TO_C,
                    rates.getId() + "." + A_TO_G,
                    rates.getId() + "." + A_TO_T,
                    rates.getId() + "." + C_TO_G,
                    rates.getId() + "." + C_TO_T,
                    rates.getId() + "." + G_TO_T});

            return new GTR(rates, freqModel);
        } else {

            Variable<Double> rateACVariable = null;
            if (xo.hasChildNamed(A_TO_C)) {
                rateACVariable = (Variable<Double>) xo.getElementFirstChild(A_TO_C);
            }
            Variable<Double> rateAGVariable = null;
            if (xo.hasChildNamed(A_TO_G)) {
                rateAGVariable = (Variable<Double>) xo.getElementFirstChild(A_TO_G);
            }
            Variable<Double> rateATVariable = null;
            if (xo.hasChildNamed(A_TO_T)) {
                rateATVariable = (Variable<Double>) xo.getElementFirstChild(A_TO_T);
            }
            Variable<Double> rateCGVariable = null;
            if (xo.hasChildNamed(C_TO_G)) {
                rateCGVariable = (Variable<Double>) xo.getElementFirstChild(C_TO_G);
            }
            Variable<Double> rateCTVariable = null;
            if (xo.hasChildNamed(C_TO_T)) {
                rateCTVariable = (Variable<Double>) xo.getElementFirstChild(C_TO_T);
            }
            Variable<Double> rateGTVariable = null;
            if (xo.hasChildNamed(G_TO_T)) {
                rateGTVariable = (Variable<Double>) xo.getElementFirstChild(G_TO_T);
            }
            int countNull = 0;
            if (rateACVariable == null) countNull++;
            if (rateAGVariable == null) countNull++;
            if (rateATVariable == null) countNull++;
            if (rateCGVariable == null) countNull++;
            if (rateCTVariable == null) countNull++;
            if (rateGTVariable == null) countNull++;

            if (countNull != 1) {
                throw new XMLParseException("Only five parameters may be specified in GTR, leave exactly one out, the others will be specifed relative to the one left out.");
            }
            return new GTR(rateACVariable, rateAGVariable, rateATVariable, rateCGVariable, rateCTVariable, rateGTVariable, freqModel);
        }
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
        return GTR.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(FREQUENCIES,
                    new XMLSyntaxRule[]{new ElementRule(FrequencyModel.class)}),
            new ElementRule(RATES,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, true),
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

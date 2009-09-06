/*
 * GTR.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
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

package dr.evomodel.substmodel;

import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.xml.*;

/**
 * General Time Reversible model of nucleotide evolution
 * This is really just a place-holder because all the implementation
 * already exists in NucleotideModel and GeneralModel, its base classes.
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: GTR.java,v 1.19 2005/05/24 20:25:58 rambaut Exp $
 */
public class GTR extends AbstractNucleotideModel {
    public static final String GTR_MODEL = "gtrModel";

    public static final String A_TO_C = "rateAC";
    public static final String A_TO_G = "rateAG";
    public static final String A_TO_T = "rateAT";
    public static final String C_TO_G = "rateCG";
    public static final String C_TO_T = "rateCT";
    public static final String G_TO_T = "rateGT";

    public static final String FREQUENCIES = "frequencies";

    private Variable<Double> rateACValue = null;
    private Variable<Double> rateAGValue = null;
    private Variable<Double> rateATValue = null;
    private Variable<Double> rateCGValue = null;
    private Variable<Double> rateCTValue = null;
    private Variable<Double> rateGTValue = null;

    /**
     * @param rateACValue rate of A<->C substitutions
     * @param rateAGValue rate of A<->G substitutions
     * @param rateATValue rate of A<->T substitutions
     * @param rateCGValue rate of C<->G substitutions
     * @param rateCTValue rate of C<->T substitutions
     * @param rateGTValue rate of G<->T substitutions
     * @param freqModel       frequencies
     */
    public GTR(
            Variable rateACValue,
            Variable rateAGValue,
            Variable rateATValue,
            Variable rateCGValue,
            Variable rateCTValue,
            Variable rateGTValue,
            FrequencyModel freqModel) {

        super(GTR_MODEL, freqModel);

        if (rateACValue != null) {
            addVariable(rateACValue);
            rateACValue.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
            this.rateACValue = rateACValue;
        }

        if (rateAGValue != null) {
            addVariable(rateAGValue);
            rateAGValue.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
            this.rateAGValue = rateAGValue;
        }

        if (rateATValue != null) {
            addVariable(rateATValue);
            rateATValue.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
            this.rateATValue = rateATValue;
        }

        if (rateCGValue != null) {
            addVariable(rateCGValue);
            rateCGValue.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
            this.rateCGValue = rateCGValue;
        }

        if (rateCTValue != null) {
            addVariable(rateCTValue);
            rateCTValue.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
            this.rateCTValue = rateCTValue;
        }

        if (rateGTValue != null) {
            addVariable(rateGTValue);
            rateGTValue.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
            this.rateGTValue = rateGTValue;
        }

    }

    public void setAbsoluteRates(double[] rates, int relativeTo) {
        for (int i = 0; i < relativeRates.length; i++) {
            relativeRates[i] = rates[i] / rates[relativeTo];
        }
        updateMatrix = true;
        fireModelChanged();
    }

    public void setRelativeRates(double[] rates) {
        System.arraycopy(rates, 0, relativeRates, 0, relativeRates.length);
        updateMatrix = true;
        fireModelChanged();
    }

    protected void setupRelativeRates() {

        if (rateACValue != null) {
            relativeRates[0] = rateACValue.getValue(0);
        }
        if (rateAGValue != null) {
            relativeRates[1] = rateAGValue.getValue(0);
        }
        if (rateATValue != null) {
            relativeRates[2] = rateATValue.getValue(0);
        }
        if (rateCGValue != null) {
            relativeRates[3] = rateCGValue.getValue(0);
        }
        if (rateCTValue != null) {
            relativeRates[4] = rateCTValue.getValue(0);
        }
        if (rateGTValue != null) {
            relativeRates[5] = rateGTValue.getValue(0);
        }
    }

    // **************************************************************
    // XHTMLable IMPLEMENTATION
    // **************************************************************

    public String toXHTML() {
        StringBuffer buffer = new StringBuffer();

        buffer.append("<em>GTR Model</em> Instantaneous Rate Matrix = <table><tr><td></td><td>A</td><td>C</td><td>G</td><td>T</td></tr>");
        buffer.append("<tr><td>A</td><td></td><td>");
        buffer.append(relativeRates[0]);
        buffer.append("</td><td>");
        buffer.append(relativeRates[1]);
        buffer.append("</td><td>");
        buffer.append(relativeRates[2]);
        buffer.append("</td></tr>");

        buffer.append("<tr><td>C</td><td></td><td></td><td>");
        buffer.append(relativeRates[3]);
        buffer.append("</td><td>");
        buffer.append(relativeRates[4]);
        buffer.append("</td></tr>");

        buffer.append("<tr><td>G</td><td></td><td></td><td></td><td>");
        buffer.append(relativeRates[5]);
        buffer.append("</td></tr>");

        buffer.append("<tr><td>G</td><td></td><td></td><td></td><td></td></tr></table>");

        return buffer.toString();
    }

    /**
     * Parses an element from an DOM document into a DemographicModel. Recognises
     * ConstantPopulation and ExponentialGrowth.
     */
    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

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
    };

}

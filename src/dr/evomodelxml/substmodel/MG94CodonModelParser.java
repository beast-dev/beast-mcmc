/*
 * MG94CodonModelParser.java
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
import dr.evomodel.substmodel.codon.MG94CodonModel;
import dr.evomodel.substmodel.codon.MG94HKYCodonModel;
import dr.evolution.datatype.Codons;
import dr.evolution.datatype.GeneticCode;
import dr.evomodel.substmodel.nucleotide.GTR;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * @author Marc A. Suchard
 * @author Guy Baele
 * @author Philippe Lemey
 */
public class MG94CodonModelParser extends AbstractXMLObjectParser {

    public static final String MUSE_CODON_MODEL = "museGautCodonModel";
    public static final String ALPHA = "alpha";
    public static final String BETA = "beta";
    public static final String KAPPA = GY94CodonModelParser.KAPPA;
    public static final String NORMALIZED = ComplexSubstitutionModelParser.NORMALIZED;

    //for GTR extension
    public static final String GTR_MODEL = GTRParser.GTR_MODEL;
    public static final String A_TO_C = GTRParser.A_TO_C;
    public static final String A_TO_G = GTRParser.A_TO_G;
    public static final String A_TO_T = GTRParser.A_TO_T;
    public static final String C_TO_G = GTRParser.C_TO_G;
    public static final String C_TO_T = GTRParser.C_TO_T;
    public static final String G_TO_T = GTRParser.G_TO_T;


    public String getParserName() {
        return MUSE_CODON_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Codons codons = Codons.UNIVERSAL;
        if (xo.hasAttribute(GeneticCode.GENETIC_CODE)) {
            String codeStr = xo.getStringAttribute(GeneticCode.GENETIC_CODE);
            codons = Codons.findByName(codeStr);
        }

        Parameter alphaParam = (Parameter) xo.getElementFirstChild(ALPHA);
        Parameter betaParam = (Parameter) xo.getElementFirstChild(BETA);
        FrequencyModel freqModel = (FrequencyModel) xo.getChild(FrequencyModel.class);

        MG94CodonModel codonModel;
        if (xo.hasChildNamed(GTR_MODEL)) {
            //TODO: change this into constructing a MG94CodonModel (needs to be written), which is started underneath
            codonModel = new MG94CodonModel(codons, alphaParam, betaParam, freqModel);

            Parameter rateACValue = null;
            if (xo.hasChildNamed(A_TO_C)) {
                rateACValue = (Parameter) xo.getElementFirstChild(A_TO_C);
            }
            Parameter rateAGValue = null;
            if (xo.hasChildNamed(A_TO_G)) {
                rateAGValue = (Parameter) xo.getElementFirstChild(A_TO_G);
            }
            Parameter rateATValue = null;
            if (xo.hasChildNamed(A_TO_T)) {
                rateATValue = (Parameter) xo.getElementFirstChild(A_TO_T);
            }
            Parameter rateCGValue = null;
            if (xo.hasChildNamed(C_TO_G)) {
                rateCGValue = (Parameter) xo.getElementFirstChild(C_TO_G);
            }
            Parameter rateCTValue = null;
            if (xo.hasChildNamed(C_TO_T)) {
                rateCTValue = (Parameter) xo.getElementFirstChild(C_TO_T);
            }
            Parameter rateGTValue = null;
            if (xo.hasChildNamed(G_TO_T)) {
                rateGTValue = (Parameter) xo.getElementFirstChild(G_TO_T);
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
            //TODO: turn this on below on and delete kappa.
            //return new new MG94GTRCodonModel(codons, alphaParam, betaParam, rateACValue, rateAGValue, rateATValue, rateCGValue, rateCTValue, rateGTValue, freqModel);

            if (xo.hasAttribute(KAPPA)){
                System.err.print("using GTR rates -- overrides KAPPA");
            }

        }  else if (xo.hasChildNamed(KAPPA)) {
            Parameter kappaParam = (Parameter)xo.getElementFirstChild(KAPPA);
            codonModel = new MG94HKYCodonModel(codons, alphaParam, betaParam, kappaParam, freqModel);
//            System.err.println("setting up MG94HKYCodonModel");
        }  else {
            //resort to standard MG94 without nucleotide rate bias
            codonModel = new MG94CodonModel(codons, alphaParam, betaParam, freqModel);
        }


        if (!xo.getAttribute(NORMALIZED, true)) {
//            codonModel.setNormalization(false);
//            Logger.getLogger("dr.app.beagle.evomodel").info("MG94CodonModel normalization: false");
        }

        return codonModel;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element represents the Yang model of codon evolution.";
    }

    public Class getReturnType() {
        return MG94CodonModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new StringAttributeRule(GeneticCode.GENETIC_CODE,
                    "The genetic code to use",
                    GeneticCode.GENETIC_CODE_NAMES, true),
            new ElementRule(ALPHA,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(BETA,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(KAPPA,
                    new XMLSyntaxRule[] { new ElementRule(Parameter.class) }, true),
            new ElementRule(A_TO_C,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, true),
            new ElementRule(A_TO_G,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, true),
            new ElementRule(A_TO_T,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, true),
            new ElementRule(C_TO_G,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, true),
            new ElementRule(C_TO_T,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, true),
            new ElementRule(G_TO_T,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, true),

            new ElementRule(FrequencyModel.class),
            AttributeRule.newBooleanRule(NORMALIZED, true),
    };

}
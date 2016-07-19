/*
 * YangCodonModelParser.java
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

import dr.evolution.datatype.Codons;
import dr.evolution.datatype.GeneticCode;
import dr.oldevomodel.substmodel.FrequencyModel;
import dr.oldevomodel.substmodel.YangCodonModel;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * Reads a the model from an XMLObject.
 */
public class YangCodonModelParser extends AbstractXMLObjectParser {
    public static final String YANG_CODON_MODEL = "yangCodonModel";
    public static final String OMEGA = "omega";
    public static final String KAPPA = "kappa";


    public String getParserName() { return YANG_CODON_MODEL; }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Codons codons = Codons.UNIVERSAL;
        if (xo.hasAttribute(GeneticCode.GENETIC_CODE)) {
            String codeStr = xo.getStringAttribute(GeneticCode.GENETIC_CODE);
            if (codeStr.equals(GeneticCode.UNIVERSAL.getName())) {
                codons = Codons.UNIVERSAL;
            } else if (codeStr.equals(GeneticCode.VERTEBRATE_MT.getName())) {
                codons = Codons.VERTEBRATE_MT;
            } else if (codeStr.equals(GeneticCode.YEAST.getName())) {
                codons = Codons.YEAST;
            } else if (codeStr.equals(GeneticCode.MOLD_PROTOZOAN_MT.getName())) {
                codons = Codons.MOLD_PROTOZOAN_MT;
            } else if (codeStr.equals(GeneticCode.INVERTEBRATE_MT.getName())) {
                codons = Codons.INVERTEBRATE_MT;
            } else if (codeStr.equals(GeneticCode.CILIATE.getName())) {
                codons = Codons.CILIATE;
            } else if (codeStr.equals(GeneticCode.ECHINODERM_MT.getName())) {
                codons = Codons.ECHINODERM_MT;
            } else if (codeStr.equals(GeneticCode.EUPLOTID_NUC.getName())) {
                codons = Codons.EUPLOTID_NUC;
            } else if (codeStr.equals(GeneticCode.BACTERIAL.getName())) {
                codons = Codons.BACTERIAL;
            } else if (codeStr.equals(GeneticCode.ALT_YEAST.getName())) {
                codons = Codons.ALT_YEAST;
            } else if (codeStr.equals(GeneticCode.ASCIDIAN_MT.getName())) {
                codons = Codons.ASCIDIAN_MT;
            } else if (codeStr.equals(GeneticCode.FLATWORM_MT.getName())) {
                codons = Codons.FLATWORM_MT;
            } else if (codeStr.equals(GeneticCode.BLEPHARISMA_NUC.getName())) {
                codons = Codons.BLEPHARISMA_NUC;
            } else if (codeStr.equals(GeneticCode.NO_STOPS.getName())) {
                codons = Codons.NO_STOPS;
            }
        }

        Parameter omegaParam = (Parameter)xo.getElementFirstChild(OMEGA);
        Parameter kappaParam = (Parameter)xo.getElementFirstChild(KAPPA);
        FrequencyModel freqModel = (FrequencyModel)xo.getChild(FrequencyModel.class);

//            codonModel.printRateMap();

        return new YangCodonModel(codons, omegaParam, kappaParam, freqModel);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element represents the Yang model of codon evolution.";
    }

    public Class getReturnType() { return YangCodonModel.class; }

    public XMLSyntaxRule[] getSyntaxRules() { return rules; }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
        new StringAttributeRule(GeneticCode.GENETIC_CODE,
            "The genetic code to use",
            new String[] {
                GeneticCode.UNIVERSAL.getName(),
                GeneticCode.VERTEBRATE_MT.getName(),
                GeneticCode.YEAST.getName(),
                GeneticCode.MOLD_PROTOZOAN_MT.getName(),
                GeneticCode.INVERTEBRATE_MT.getName(),
                GeneticCode.CILIATE.getName(),
                GeneticCode.ECHINODERM_MT.getName(),
                GeneticCode.EUPLOTID_NUC.getName(),
                GeneticCode.BACTERIAL.getName(),
                GeneticCode.ALT_YEAST.getName(),
                GeneticCode.ASCIDIAN_MT.getName(),
                GeneticCode.FLATWORM_MT.getName(),
                GeneticCode.BLEPHARISMA_NUC.getName(),
                GeneticCode.NO_STOPS.getName()}, true),
        new ElementRule(OMEGA,
            new XMLSyntaxRule[] { new ElementRule(Parameter.class) }),
        new ElementRule(KAPPA,
            new XMLSyntaxRule[] { new ElementRule(Parameter.class) }),
        new ElementRule(FrequencyModel.class)
    };
}

/*
 * EmpiricalAminoAcidModelParser.java
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

import dr.evomodel.substmodel.EmpiricalRateMatrix;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.aminoacid.*;
import dr.xml.*;

/**
 * BEAGLE parser for empirical amino acid models
 */
public class EmpiricalAminoAcidModelParser extends AbstractXMLObjectParser {

    public static final String EMPIRICAL_AMINO_ACID_MODEL = "aminoAcidModel";
    public static final String FREQUENCIES = "frequencies";
    public static final String TYPE = "type";

    public String getParserName() {
        return EMPIRICAL_AMINO_ACID_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        FrequencyModel freqModel = null;

        if (xo.hasAttribute(FREQUENCIES)) {
            XMLObject cxo = xo.getChild(FREQUENCIES);
            freqModel = (FrequencyModel) cxo.getChild(FrequencyModel.class);
        }

        EmpiricalRateMatrix rateMatrix = null;

        String type = xo.getStringAttribute(TYPE);

        if (type.equals(AminoAcidModelType.BLOSUM_62.getXMLName())) {
            rateMatrix = Blosum62.INSTANCE;
        } else if (type.equals(AminoAcidModelType.DAYHOFF.getXMLName())) {
            rateMatrix = Dayhoff.INSTANCE;
        } else if (type.equals(AminoAcidModelType.JTT.getXMLName())) {
            rateMatrix = JTT.INSTANCE;
        } else if (type.equals(AminoAcidModelType.MT_REV_24.getXMLName())) {
            rateMatrix = MTREV.INSTANCE;
        } else if (type.equals(AminoAcidModelType.CP_REV_45.getXMLName())) {
            rateMatrix = CPREV.INSTANCE;
        } else if (type.equals(AminoAcidModelType.WAG.getXMLName())) {
            rateMatrix = WAG.INSTANCE;
        } else if (type.equals(AminoAcidModelType.LG.getXMLName())) {
            rateMatrix = LG.INSTANCE;
        } else if (type.equals(AminoAcidModelType.FLU.getXMLName())) {
            rateMatrix = FLU.INSTANCE;
        } else {
            throw new XMLParseException("Unrecognized empirical amino acid model: " + type);
        }

        return new EmpiricalAminoAcidModel(rateMatrix, freqModel);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new StringAttributeRule(TYPE, "The type of empirical amino-acid rate matrix", AminoAcidModelType.xmlNames(), false),
            new ElementRule(FREQUENCIES, FrequencyModel.class, "If the frequencies are omitted than the empirical frequencies associated with the selected model are used.", true)
    };

    public String getParserDescription() {
        return "An empirical amino acid substitution model.";
    }

    public Class getReturnType() {
        return EmpiricalAminoAcidModel.class;
    }
}
/*
 * ConvertAlignmentParser.java
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

package dr.evoxml;

import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.ConvertAlignment;
import dr.evolution.datatype.*;
import dr.xml.*;

import java.util.logging.Logger;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 *
 * @version $Id: ConvertAlignmentParser.java,v 1.3 2005/07/11 14:06:25 rambaut Exp $
 */
public class ConvertAlignmentParser extends AbstractXMLObjectParser {

    public final static String CONVERT = "convert";

    public String getParserName() { return CONVERT; }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Alignment alignment = (Alignment)xo.getChild(Alignment.class);

        DataType dataType = null;

        String dataTypeStr = xo.getStringAttribute(DataType.DATA_TYPE);

        if (dataTypeStr.equals(Nucleotides.DESCRIPTION)) {
            dataType = Nucleotides.INSTANCE;
        } else if (dataTypeStr.equals(AminoAcids.DESCRIPTION)) {
            dataType = AminoAcids.INSTANCE;
        } else if (dataTypeStr.equals(Codons.DESCRIPTION)) {
            dataType = Codons.UNIVERSAL;
        } else if (dataTypeStr.equals(TwoStates.DESCRIPTION)) {
            dataType = TwoStates.INSTANCE;
        }

        GeneticCode geneticCode = GeneticCode.UNIVERSAL;
        if (xo.hasAttribute(GeneticCode.GENETIC_CODE)) {
            String codeStr = xo.getStringAttribute(GeneticCode.GENETIC_CODE);
            if (codeStr.equals(GeneticCode.UNIVERSAL.getName())) {
                geneticCode = GeneticCode.UNIVERSAL;
            } else if (codeStr.equals(GeneticCode.VERTEBRATE_MT.getName())) {
                geneticCode = GeneticCode.VERTEBRATE_MT;
            } else if (codeStr.equals(GeneticCode.YEAST.getName())) {
                geneticCode = GeneticCode.YEAST;
            } else if (codeStr.equals(GeneticCode.MOLD_PROTOZOAN_MT.getName())) {
                geneticCode = GeneticCode.MOLD_PROTOZOAN_MT;
            } else if (codeStr.equals(GeneticCode.INVERTEBRATE_MT.getName())) {
                geneticCode = GeneticCode.INVERTEBRATE_MT;
            } else if (codeStr.equals(GeneticCode.CILIATE.getName())) {
                geneticCode = GeneticCode.CILIATE;
            } else if (codeStr.equals(GeneticCode.ECHINODERM_MT.getName())) {
                geneticCode = GeneticCode.ECHINODERM_MT;
            } else if (codeStr.equals(GeneticCode.EUPLOTID_NUC.getName())) {
                geneticCode = GeneticCode.EUPLOTID_NUC;
            } else if (codeStr.equals(GeneticCode.BACTERIAL.getName())) {
                geneticCode = GeneticCode.BACTERIAL;
            } else if (codeStr.equals(GeneticCode.ALT_YEAST.getName())) {
                geneticCode = GeneticCode.ALT_YEAST;
            } else if (codeStr.equals(GeneticCode.ASCIDIAN_MT.getName())) {
                geneticCode = GeneticCode.ASCIDIAN_MT;
            } else if (codeStr.equals(GeneticCode.FLATWORM_MT.getName())) {
                geneticCode = GeneticCode.FLATWORM_MT;
            } else if (codeStr.equals(GeneticCode.BLEPHARISMA_NUC.getName())) {
                geneticCode = GeneticCode.BLEPHARISMA_NUC;
            }
        }

        ConvertAlignment convert = new ConvertAlignment(dataType, geneticCode, alignment);
	    Logger.getLogger("dr.evoxml").info("Converted alignment, '" + xo.getId() + "', from " +
	            alignment.getDataType().getDescription() + " to " + dataType.getDescription());


        return convert;
    }

    public String getParserDescription() {
        return "Converts an alignment to the given data type.";
    }

    public Class getReturnType() { return Alignment.class; }

    public XMLSyntaxRule[] getSyntaxRules() { return rules; }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
        new ElementRule(Alignment.class),
        new StringAttributeRule(DataType.DATA_TYPE,
            "The type of sequence data",
            new String[] {Nucleotides.DESCRIPTION, AminoAcids.DESCRIPTION, Codons.DESCRIPTION, TwoStates.DESCRIPTION},
            false )
    };
}

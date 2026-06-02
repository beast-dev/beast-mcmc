/*
 * ConvertAlignmentParser.java
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

package dr.evoxml;

import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.ConvertAlignment;
import dr.evolution.datatype.*;
import dr.evoxml.util.DataTypeUtils;
import dr.xml.*;

import java.util.logging.Logger;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 *
 */
public class ConvertAlignmentParser extends AbstractXMLObjectParser {

    public final static String CONVERT = "convert";

    public String getParserName() { return CONVERT; }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Alignment alignment = (Alignment)xo.getChild(Alignment.class);

        DataType dataType = (DataType)xo.getChild(DataType.class);
        if (dataType == null) {
            dataType = DataTypeUtils.getDataType(xo);
        }

        GeneticCode geneticCode = GeneticCode.UNIVERSAL;
        if (dataType instanceof Codons) {
            geneticCode = ((Codons)dataType).getGeneticCode();
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
            new XORRule(new XMLSyntaxRule[]{
                    new ElementRule(DataType.class),
                    new StringAttributeRule(DataType.DATA_TYPE,
                            "The type of sequence data",
                            new String[]{Nucleotides.DESCRIPTION, AminoAcids.DESCRIPTION, Codons.DESCRIPTION, TwoStates.DESCRIPTION,
                                    HiddenCodons.DESCRIPTION + "2", HiddenCodons.DESCRIPTION + "3", HiddenCodons.DESCRIPTION + "4",
                                    HiddenCodons.DESCRIPTION + "5"}, false)
            })
    };
}

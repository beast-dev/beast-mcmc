/*
 * HypermutantAlignmentParser.java
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

package dr.evoxml;

import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.HypermutantAlignment;
import dr.evolution.datatype.*;
import dr.xml.*;

import java.util.logging.Logger;

/**
 * This converts 'A's that are at specific APOBEC targeted contexts into an A/G ambiguity code for
 * later recognition by the HypermutantErrorModel.
 *
 * @author Andrew Rambaut
 *
 * @version $Id$
 */
public class HypermutantAlignmentParser extends AbstractXMLObjectParser {
    public final static String HYPERMUTANT_ALIGNMENT = "hypermutantAlignment";
    public final static String CONTEXT_TYPE = "type";

    public String getParserName() { return HYPERMUTANT_ALIGNMENT; }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Alignment alignment = (Alignment)xo.getChild(Alignment.class);

        if (alignment.getDataType().getType() != DataType.NUCLEOTIDES) {
            throw new XMLParseException("HypermutantAlignment can only convert nucleotide alignments");
        }

        String typeName = xo.getStringAttribute(CONTEXT_TYPE);
        HypermutantAlignment.APOBECType type = null;
        try {
            type = HypermutantAlignment.APOBECType.valueOf(typeName.toUpperCase());
        } catch(IllegalArgumentException iae) {
            throw new XMLParseException("Unrecognised hypermutation type: " + typeName);
        }

        HypermutantAlignment convert = new HypermutantAlignment(type, alignment);
        int mutatedCount = convert.getMutatedContextCount();
        int totalCount = mutatedCount + convert.getUnmutatedContextCount();

        Logger.getLogger("dr.evoxml").info("Converted alignment, '" + xo.getId() + "' to a hypermutant alignment targeting " + type.toString() + " contexts.\r" +
                "\tPotentially mutated contexts: " + mutatedCount + " out of a total of " + totalCount + " contexts");

        return convert;
    }

    public String getParserDescription() {
        return "Converts an alignment so that 'A's at specific APOBEC targeted contexts are set to an A/G ambiguity code.";
    }

    public Class getReturnType() { return Alignment.class; }

    public XMLSyntaxRule[] getSyntaxRules() { return rules; }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
            new ElementRule(Alignment.class),
            new StringAttributeRule(CONTEXT_TYPE,
                    "The type of APOBEC molecule being modelled",
                    new String[] {HypermutantAlignment.APOBECType.ALL.toString(), HypermutantAlignment.APOBECType.BOTH.toString(), HypermutantAlignment.APOBECType.HA3G.toString(), HypermutantAlignment.APOBECType.HA3F.toString()},
                    false )
    };
}
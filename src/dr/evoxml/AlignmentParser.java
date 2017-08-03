/*
 * AlignmentParser.java
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
import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.sequence.Sequence;
import dr.evolution.sequence.UncertainSequence;
import dr.evoxml.util.DataTypeUtils;
import dr.xml.*;

import java.util.logging.Logger;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @version $Id: AlignmentParser.java,v 1.3 2005/07/11 14:06:25 rambaut Exp $
 */
public class AlignmentParser extends AbstractXMLObjectParser {

    public static final String ALIGNMENT = "alignment";

    public String getParserName() {
        return ALIGNMENT;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        final SimpleAlignment alignment = new SimpleAlignment();

        final DataType dataType = DataTypeUtils.getDataType(xo);

        if (dataType == null) {
            throw new XMLParseException("dataType attribute expected for alignment element");
        }

        alignment.setDataType(dataType);

        for (int i = 0; i < xo.getChildCount(); i++) {

            final Object child = xo.getChild(i);
            if (child instanceof UncertainSequence) {
                alignment.addSequence((UncertainSequence) child);
            } else if (child instanceof Sequence) {
                alignment.addSequence((Sequence) child);
            } else if (child instanceof DataType) {
                // already dealt with
            } else {
                throw new XMLParseException("Unknown child element found in alignment");
            }
        }

        final Logger logger = Logger.getLogger("dr.evoxml");
        logger.info("\nRead alignment" + (xo.hasAttribute(XMLParser.ID) ? ": " + xo.getId() : "") +
                "\n  Sequences = " + alignment.getSequenceCount() +
                "\n      Sites = " + alignment.getSiteCount() +
                "\n   Datatype = " + alignment.getDataType().getDescription());

        return alignment;
    }

    public String getParserDescription() {
        return "This element represents an alignment of molecular sequences.";
    }

    public Class getReturnType() {
        return Alignment.class;
    }

    public String getExample() {

        return
                "<!-- An alignment of three short DNA sequences -->\n" +
                        "<alignment missing=\"-?\" dataType=\"" + Nucleotides.DESCRIPTION + "\">\n" +
                        "  <sequence>\n" +
                        "    <taxon idref=\"taxon1\"/>\n" +
                        "    ACGACTAGCATCGAGCTTCG--GATAGCAGGC\n" +
                        "  </sequence>\n" +
                        "  <sequence>\n" +
                        "    <taxon idref=\"taxon2\"/>\n" +
                        "    ACGACTAGCATCGAGCTTCGG-GATAGCATGC\n" +
                        "  </sequence>\n" +
                        "  <sequence>\n" +
                        "    <taxon idref=\"taxon3\"/>\n" +
                        "    ACG?CTAGAATCGAGCTTCGAGGATAGCATGC\n" +
                        "  </sequence>\n" +
                        "</alignment>\n";
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new XORRule(
                    new StringAttributeRule(
                            DataType.DATA_TYPE,
                            "The data type",
                            DataType.getRegisteredDataTypeNames(), false),
                    new ElementRule(DataType.class)
            ),
            new ElementRule(Sequence.class, 1, Integer.MAX_VALUE)
    };
}

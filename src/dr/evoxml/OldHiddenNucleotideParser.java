/*
 * HiddenNucleotideParser.java
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

import dr.evolution.datatype.OldHiddenNucleotides;
import dr.xml.*;

/**
 * @author Alexei Drummond
 */
@Deprecated
public class OldHiddenNucleotideParser extends AbstractXMLObjectParser {

    public static final String HIDDEN_NUCLEOTIDES = "hiddenNucleotides";
    public static final String HIDDEN_CLASS_COUNT = "classCount";

    public String getParserName() {
        return HIDDEN_NUCLEOTIDES;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        int hiddenClassCount = xo.getIntegerAttribute(HIDDEN_CLASS_COUNT);
        return new OldHiddenNucleotides(hiddenClassCount);
    }

    public String getParserDescription() {
        return "A nucleotide data type that allows hidden substitution classes";
    }

    public Class getReturnType() {
        return OldHiddenNucleotides.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{AttributeRule.newIntegerRule(HIDDEN_CLASS_COUNT)};
    }
}
/*
 * ExpandAlignmentStateSpaceParser.java
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
import dr.evolution.alignment.ConvertAlignment;
import dr.evolution.alignment.ExpandAlignmentStateSpace;
import dr.evolution.alignment.PatternList;
import dr.evolution.datatype.*;
import dr.evolution.util.Taxa;
import dr.evolution.util.TaxonList;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.treedatalikelihood.discrete.SequenceDistanceStatistic;
import dr.evoxml.util.DataTypeUtils;
import dr.xml.*;

import java.util.logging.Logger;

public class ExpandAlignmentStateSpaceParser extends AbstractXMLObjectParser {

    public final static String EXPAND = "expandStateSpace";
    public final static String TRAIT = "trait";
    public final static String ALIGNMENT = "alignment";

    public String getParserName() { return EXPAND; }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        DataType dataType = DataTypeUtils.getDataType(xo);

        Alignment alignment = (Alignment)xo.getElementFirstChild(ALIGNMENT);

        PatternList trait = (PatternList)xo.getElementFirstChild(TRAIT);

        ExpandAlignmentStateSpace expand = new ExpandAlignmentStateSpace(dataType, alignment, trait);
        Logger.getLogger("dr.evoxml").info("Converted alignment, '" + xo.getId() + "', from " +
                alignment.getDataType().getDescription() + " to " + dataType.getDescription());


        return expand;
    }

    public String getParserDescription() {
        return "Converts an alignment to the given data type.";
    }

    public Class getReturnType() { return Alignment.class; }

    public XMLSyntaxRule[] getSyntaxRules() { return rules; }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
            new ElementRule(ALIGNMENT,Alignment.class),
            new ElementRule(TRAIT,PatternList.class)
    };
}

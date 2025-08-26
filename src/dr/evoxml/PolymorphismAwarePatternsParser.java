/*
 * PolymorphismAwarePatternsParser.java
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

import dr.evolution.alignment.*;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.PolymorphismAwareDataType;
import dr.evolution.sequence.Sequence;
import dr.evolution.sequence.UncertainSequence;
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.xml.*;

import java.util.List;

import static dr.evoxml.SitePatternsParser.STRIP;
import static dr.evoxml.SitePatternsParser.FROM;
import static dr.evoxml.SitePatternsParser.TO;
import static dr.evoxml.SitePatternsParser.EVERY;

/**
 * @author Xiang Ji
 * @author Nicola De Maio
 * @author Ben Redelings
 * @author Marc A. Suchard
 */
public class PolymorphismAwarePatternsParser extends AbstractXMLObjectParser {

    public static final String NAME = "polymorphismAwarePatterns";
    public static final String VIRTUAL_POP_SIZE = "virtualPopSize";


    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        SimpleAlignment alignment = (SimpleAlignment) xo.getChild(SimpleAlignment.class);
        DataType baseDataType = alignment.getDataType();
        int virtualPopSize = xo.getIntegerAttribute(VIRTUAL_POP_SIZE);
        Taxa taxa = (Taxa) xo.getChild(Taxa.class);

        PolymorphismAwareDataType dataType = new PolymorphismAwareDataType(baseDataType, virtualPopSize);
        DataType.registerDataType(dataType.getDataTypeDescription(baseDataType, virtualPopSize), dataType);

        int from = 0;
        int to = -1;
        int every = xo.getAttribute(EVERY, 1);

        boolean strip = xo.getAttribute(STRIP, true);

        if (xo.hasAttribute(FROM)) {
            from = xo.getIntegerAttribute(FROM) - 1;

            if (from < 0)
                throw new XMLParseException("illegal 'from' attribute in patterns element");
        }

        if (xo.hasAttribute(TO)) {
            to = xo.getIntegerAttribute(TO) - 1;
            if (to < 0 || to < from)
                throw new XMLParseException("illegal 'to' attribute in patterns element");
        }

        SimpleSiteList patterns = new SimpleSiteList(dataType, taxa);

        int[][] polymorphismAwarePattern = new int[alignment.getSiteCount()][];
        for (int i = 0; i < alignment.getSiteCount(); i++) {
            polymorphismAwarePattern[i] = new int[taxa.getTaxonCount()];
        }

        for (int taxonIndex = 0; taxonIndex < taxa.getTaxonCount(); taxonIndex++) {
            UncertainSequence sequence = (UncertainSequence) alignment.getSequence(taxonIndex);
            for (int siteIndex = 0; siteIndex < alignment.getSiteCount(); siteIndex++) {
                UncertainSequence.UncertainCharacterList characters = sequence.getUncertainCharacterList(siteIndex);
                polymorphismAwarePattern[siteIndex][taxonIndex] = dataType.getState(characters);
            }
        }

        for (int i = 0; i < alignment.getSiteCount(); i++) {
            patterns.addPattern(polymorphismAwarePattern[i]);
        }

        return patterns;
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() { return rules; }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
            new ElementRule(SimpleAlignment.class),
            AttributeRule.newIntegerRule(VIRTUAL_POP_SIZE),
            new ElementRule(Taxa.class, "The taxon set")
    };

    @Override
    public String getParserDescription() {
        return "A PoMo site pattern";
    }

    @Override
    public Class getReturnType() {
        return PatternList.class;
    }

    @Override
    public String getParserName() {
        return NAME;
    }
}

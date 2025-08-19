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

import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.alignment.SitePatterns;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.PolymorphismAwareDataType;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

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
        TaxonList taxa = (TaxonList) xo.getChild(TaxonList.class);

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

        SitePatterns patterns = new SitePatterns(alignment, from, to, every, strip);

        double[][] polymorphismAwarePattern = new double[taxa.getTaxonCount()][];

        for (int i = 0; i < taxa.getTaxonCount(); i++) {
            Taxon taxon = taxa.getTaxon(i);
        }


        return null;
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[0];
    }

    @Override
    public String getParserDescription() {
        return "";
    }

    @Override
    public Class getReturnType() {
        return null;
    }

    @Override
    public String getParserName() {
        return NAME;
    }
}

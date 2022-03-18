/*
 * PairedParalogSitePatternsParser.java
 *
 * Copyright (c) 2002-2017 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

import dr.evolution.alignment.PairedParalogSitePatterns;
import dr.evolution.alignment.SitePatterns;
import dr.evolution.util.Taxa;
import dr.xml.*;

/**
 * @author Xiang Ji
 * @author Jeff Thorne
 * @author Marc A. Suchard
 */
public class PairedParalogSitePatternsParser extends AbstractXMLObjectParser {

    private final String NAME = "pairedParalogSitePatterns";
    private final String ID_SEPARATOR = "idSeparator";
    private final String PARALOGS = "paralogs";
    private final String SINGLE_COPY_SPECIES = "singleCopySpecies";
    private final String SIZE = "size";
    private final String SPECIES = "species";
    private final String ALL_SEQ = "allSeq";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        String idSeparator = (String) xo.getAttribute(ID_SEPARATOR, "__");
        String[] paralogs = xo.getStringArrayAttribute(PARALOGS);
        final int nParalogs = (int) xo.getAttribute(SIZE, 2);
        final SitePatterns siteList = (SitePatterns) xo.getChild(SitePatterns.class);
        final Taxa speciesTaxa = (Taxa) xo.getChild(SPECIES).getChild(Taxa.class);
        final String[] singleCopySpecies = xo.getStringArrayAttribute(SINGLE_COPY_SPECIES);

        if (paralogs.length != nParalogs) {
            throw new RuntimeException("Paralog list dimension mismatch input number of paralogs.");
        }

        if (nParalogs != 2) {
            throw new RuntimeException("Not yet implemented for more than two paralogs.");
        }

        return new PairedParalogSitePatterns(siteList, paralogs, idSeparator, speciesTaxa, singleCopySpecies);
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newStringRule(ID_SEPARATOR, true),
            AttributeRule.newStringArrayRule(PARALOGS),
            AttributeRule.newStringArrayRule(SINGLE_COPY_SPECIES),
            AttributeRule.newIntegerRule(SIZE, true),
            new ElementRule(SitePatterns.class),
            new ElementRule(SPECIES, Taxa.class),
    };

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return PairedParalogSitePatterns.class;
    }

    @Override
    public String getParserName() {
        return NAME;
    }
}

/*
 * RandomSubsetTaxaParser.java
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

import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.math.MathUtils;
import dr.xml.*;

import java.util.List;

/**
 * @author Marc A. Suchard
 */

public class RandomSubsetTaxaParser extends AbstractXMLObjectParser {

     public static final String RANDOM_SUBSET_TAXA = "randomSubsetTaxa";
    public static final String COUNT = "total";
    public static final String WITH_REPLACEMENT = "withReplacement";

    public String getParserName() { return RANDOM_SUBSET_TAXA; }


    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Taxa originalTaxa = new Taxa();

        for (int i = 0; i < xo.getChildCount(); i++) {
            Object child = xo.getChild(i);
            if (child instanceof Taxon) {
                Taxon taxon = (Taxon)child;
                originalTaxa.addTaxon(taxon);
            } else if (child instanceof TaxonList) {
                TaxonList taxonList1 = (TaxonList)child;
                for (int j = 0; j < taxonList1.getTaxonCount(); j++) {
                    originalTaxa.addTaxon(taxonList1.getTaxon(j));
                }
            } else {
                throwUnrecognizedElement(xo);
            }
        }

        List<Taxon> originalTaxonList = originalTaxa.asList();

        int sampleTotal = xo.getAttribute(COUNT, originalTaxonList.size());
        if (sampleTotal < 2) {
            throw new XMLParseException("Must sample atleast two taxa");
        }

        boolean withReplacement = xo.getAttribute(WITH_REPLACEMENT, false);

        Taxa sampledTaxa = new Taxa();
        for (int i = 0; i < sampleTotal; i++) {
            int choice = MathUtils.nextInt(originalTaxonList.size());
            Taxon taxonToAdd = originalTaxonList.get(choice);
            sampledTaxa.addTaxon(taxonToAdd);

            if (!withReplacement) {
                originalTaxonList.remove(choice);
            }
        }

        return sampledTaxa;
    }

    public XMLSyntaxRule[] getSyntaxRules() { return rules; }

    private final XMLSyntaxRule[] rules = {
        new OrRule(
            new ElementRule(Taxa.class, 1, Integer.MAX_VALUE),
            new ElementRule(Taxon.class, 1, Integer.MAX_VALUE)
        ),
            AttributeRule.newIntegerRule(COUNT),
            AttributeRule.newBooleanRule(WITH_REPLACEMENT),
    };

    public String getParserDescription() {
        return "Defines a set of taxon objects.";
    }

    public Class getReturnType() { return Taxa.class; }
}

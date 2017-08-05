/*
 * TaxaParser.java
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
import dr.xml.*;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 *
 * @version $Id: TaxaParser.java,v 1.2 2005/05/24 20:25:59 rambaut Exp $
 */
public class TaxaParser extends AbstractXMLObjectParser {

    public static final String TAXA = "taxa";
    public static final String EXCLUDE = "exclude";

    public String getParserName() { return TAXA; }

    public String getExample() {
        return "<!-- A list of six taxa -->\n"+
                "<taxa id=\"greatApes\">\n"+
                "	<taxon id=\"human\"/>\n"+
                "	<taxon id=\"chimp\"/>\n"+
                "	<taxon id=\"bonobo\"/>\n"+
                "	<taxon id=\"gorilla\"/>\n"+
                "	<taxon id=\"orangutan\"/>\n"+
                "	<taxon id=\"siamang\"/>\n"+
                "</taxa>\n" +
                "\n" +
                "<!-- A list of three taxa by references to above taxon objects -->\n"+
                "<taxa id=\"humanAndChimps\">\n"+
                "	<taxon idref=\"human\"/>\n"+
                "	<taxon idref=\"chimp\"/>\n"+
                "	<taxon idref=\"bonobo\"/>\n"+
                "</taxa>\n";
    }

    /** @return an instance of Node created from a DOM element */
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Taxa taxonList = new Taxa();

        for (Taxon taxon : xo.getAllChildren(Taxon.class)) {
            taxonList.addTaxon(taxon);
        }
        
        for (Taxa taxa : xo.getAllChildren(Taxa.class)) {
            taxonList.addTaxa(taxa);
        }

        for (XMLObject cxo : xo.getAllChildren(EXCLUDE)) {
            for (Taxa exclude : cxo.getAllChildren(Taxa.class)) {
                taxonList.removeTaxa(exclude);
            }
        }

        return taxonList;
    }

    public XMLSyntaxRule[] getSyntaxRules() { return rules; }

    private final XMLSyntaxRule[] rules = {
            new OrRule(
                    new ElementRule(Taxa.class, 1, Integer.MAX_VALUE),
                    new ElementRule(Taxon.class, 1, Integer.MAX_VALUE)
            ),
            new ElementRule("exclude", Taxa.class, "taxa to exclude", 0, Integer.MAX_VALUE)
    };

    public String getParserDescription() {
        return "Defines a set of taxon objects.";
    }

    public Class getReturnType() { return Taxa.class; }
}



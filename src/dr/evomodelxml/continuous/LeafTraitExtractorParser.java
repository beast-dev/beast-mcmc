/*
 * LeafTraitExtractorParser.java
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

package dr.evomodelxml.continuous;

import dr.evolution.tree.MutableTree;
import dr.evomodelxml.tree.TreeModelParser;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * @author Marc A. Suchard
 */
public class LeafTraitExtractorParser extends AbstractXMLObjectParser {

    public static final String NAME = "leafTraitParameter";
    public static final String SET_BOUNDS = "setBounds";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        MutableTree model = (MutableTree) xo.getChild(MutableTree.class);
        final CompoundParameter allTraits = (CompoundParameter) xo.getChild(CompoundParameter.class);

        String taxonString = (String) xo.getAttribute(TreeModelParser.TAXON);
        final int leafIndex = model.getTaxonIndex(taxonString);
        if (leafIndex == -1) {
            throw new XMLParseException("Unable to find taxon '" + taxonString + "' in trees.");
        }
        final Parameter leafTrait = allTraits.getParameter(leafIndex);

        boolean setBounds = xo.getAttribute(SET_BOUNDS, true);
        if (setBounds) {

            Parameter.DefaultBounds bound = new Parameter.DefaultBounds(Double.MAX_VALUE, -Double.MAX_VALUE,
                    leafTrait.getDimension());
            leafTrait.addBounds(bound);
        }

        return leafTrait;
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                AttributeRule.newStringRule(TreeModelParser.TAXON),
                AttributeRule.newBooleanRule(SET_BOUNDS, true),
                new ElementRule(MutableTree.class),
                new ElementRule(CompoundParameter.class),
        };
    }

    @Override
    public String getParserDescription() {
        return "Parses the leaf trait parameter out of the compound parameter of an integrated trait likelihood";
    }

    @Override
    public Class getReturnType() {
        return Parameter.class;
    }

    public String getParserName() {
        return NAME;
    }
}

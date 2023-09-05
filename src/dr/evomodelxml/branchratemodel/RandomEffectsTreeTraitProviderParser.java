/*
 * RandomEffectsTreeTraitProviderParser.java
 *
 * Copyright (c) 2002-2022 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodelxml.branchratemodel;

import dr.evomodel.branchratemodel.ArbitraryBranchRates;
import dr.evomodel.branchratemodel.RandomEffectsTreeTraitProvider;
import dr.xml.*;

/**
 * @author Marc A. Suchard
 */
public class RandomEffectsTreeTraitProviderParser extends AbstractXMLObjectParser {

    private static final String PARSER_NAME = "randomEffectsTrait";
    private static final String TRAIT_NAME = "name";
    private static final String TAKE_LOG = "log";

    public String getParserName() {
        return PARSER_NAME;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        final ArbitraryBranchRates rates = (ArbitraryBranchRates) xo.getChild(ArbitraryBranchRates.class);

        final String traitName = (xo.hasAttribute(TRAIT_NAME) ?
                xo.getStringAttribute(TRAIT_NAME) : "randomEffect");

        final boolean takeLog = xo.getAttribute(TAKE_LOG, false);

        return new RandomEffectsTreeTraitProvider(rates, traitName, takeLog);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element returns a branch rate model." +
                "The branch rates are specified by an attribute embedded in the nodes of the tree.";
    }

    public Class getReturnType() {
        return RandomEffectsTreeTraitProvider.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(ArbitraryBranchRates.class),
            AttributeRule.newBooleanRule(TAKE_LOG, true),
            new StringAttributeRule(TRAIT_NAME,
                    "Optional name of a rate attribute to be read with the trees", true)
    };
}

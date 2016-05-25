/*
 * ImportanceNarrowExchangeParser.java
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

package dr.evomodelxml.operators;

import dr.evolution.alignment.PatternList;
import dr.evomodel.operators.ImportanceNarrowExchange;
import dr.evomodel.tree.TreeModel;
import dr.inference.operators.MCMCOperator;
import dr.xml.*;

/**
 */
public class ImportanceNarrowExchangeParser extends AbstractXMLObjectParser {

    public static final String INS = "ImportanceNarrowExchange";
    public static final String EPSILON = "epsilon";

    public String getParserName() {
        return INS;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        final TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
        final PatternList patterns = (PatternList) xo.getChild(PatternList.class);
        final double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
        final double epsilon = xo.getAttribute(EPSILON, 0.1);

        try {
            return new ImportanceNarrowExchange(treeModel, patterns, epsilon, weight);
        } catch( Exception e ) {
            throw new XMLParseException(e.getMessage());
        }
    }

    // ************************************************************************
    // AbstractXMLObjectParser
    // implementation
    // ************************************************************************

    public String getParserDescription() {
        return "This element represents a swap operator. "
                + "This operator swaps a random subtree with its uncle.";
    }

    public Class getReturnType() {
        return ImportanceNarrowExchange.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
            AttributeRule.newDoubleRule(EPSILON, true),
            new ElementRule(PatternList.class),
            new ElementRule(TreeModel.class)
    };
}

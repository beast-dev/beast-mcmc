/*
 * ExchangeOperatorParser.java
 *
 * Copyright (C) 2002-2011 Alexei Drummond and Andrew Rambaut
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
 * BEAST is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.evomodelxml.operators;

import dr.evomodel.operators.ExchangeOperator;
import dr.evomodel.tree.TreeModel;
import dr.inference.operators.MCMCOperator;
import dr.xml.*;

/**
 */
public class ExchangeOperatorParser {

    public static final String NARROW_EXCHANGE = "narrowExchange";
    public static final String WIDE_EXCHANGE = "wideExchange";
    public static final String INTERMEDIATE_EXCHANGE = "intermediateExchange";

    public static XMLObjectParser NARROW_EXCHANGE_OPERATOR_PARSER = new AbstractXMLObjectParser() {
        public String getParserName() {
            return NARROW_EXCHANGE;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            final TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
            final double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
            return new ExchangeOperator(ExchangeOperator.NARROW, treeModel, weight);
        }

        // ************************************************************************
        // AbstractXMLObjectParser implementation
        // ************************************************************************

        public String getParserDescription() {
            return "This element represents a narrow exchange operator. "
                    + "This operator swaps a random subtree with its uncle.";
        }

        public Class getReturnType() {
            return ExchangeOperator.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
                new ElementRule(TreeModel.class)
        };
    };

    public static XMLObjectParser INTERMEDIATE_EXCHANGE_OPERATOR_PARSER = new AbstractXMLObjectParser() {
        public String getParserName() {
            return INTERMEDIATE_EXCHANGE;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            final TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
            final double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
            return new ExchangeOperator(ExchangeOperator.INTERMEDIATE, treeModel, weight);
        }

        // ************************************************************************
        // AbstractXMLObjectParser implementation
        // ************************************************************************

        public String getParserDescription() {
            return "This element represents an intermediate exchange operator. "
                    + "This operator swaps two subtree random subtrees with a bias towards nearby subtrees.";
        }

        public Class getReturnType() {
            return ExchangeOperator.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
                new ElementRule(TreeModel.class)
        };
    };


    public static XMLObjectParser WIDE_EXCHANGE_OPERATOR_PARSER = new AbstractXMLObjectParser() {
        public String getParserName() {
            return WIDE_EXCHANGE;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            final TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
            if (treeModel.getExternalNodeCount() <= 2) {
                throw new XMLParseException("Tree with less than 3 taxa");
            }
            final double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);

            return new ExchangeOperator(ExchangeOperator.WIDE, treeModel, weight);
        }

        // ************************************************************************
        // AbstractXMLObjectParser implementation
        // ************************************************************************

        public String getParserDescription() {
            return "This element represents a wide exchange operator. "
                    + "This operator swaps two random subtrees.";
        }

        public Class getReturnType() {
            return ExchangeOperator.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules;{
            rules = new XMLSyntaxRule[]{
                    AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
                    new ElementRule(TreeModel.class)
            };
        }
    };
}

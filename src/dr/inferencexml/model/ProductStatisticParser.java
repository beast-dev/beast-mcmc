/*
 * ProductStatisticParser.java
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

package dr.inferencexml.model;

import dr.inference.model.ProductStatistic;
import dr.inference.model.Statistic;
import dr.xml.*;

/**
 */
public class ProductStatisticParser extends AbstractXMLObjectParser {

    public static String PRODUCT_STATISTIC = "productStatistic";
    public static String PRODUCT = "product";
    public static String ELEMENT_WISE = "elementwise";
    public static String CONSTANT = "constant";

    public String[] getParserNames() {
        return new String[]{getParserName(), PRODUCT};
    }

    public String getParserName() {
        return PRODUCT_STATISTIC;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {


        boolean elementwise = xo.getAttribute(ELEMENT_WISE, false);

        String name = PRODUCT_STATISTIC;
        if (xo.hasAttribute(Statistic.NAME)) {
            name = xo.getAttribute(Statistic.NAME, xo.getId());
        } else if (xo.hasAttribute(XMLParser.ID)) {
            name = xo.getAttribute(XMLParser.ID, xo.getId());
        }

        double[] constants = null;
        if (xo.hasAttribute(CONSTANT)) {
            constants = xo.getDoubleArrayAttribute(CONSTANT);
        }

        if (constants != null && constants.length != 1 && elementwise) {
            throw new XMLParseException("The constant given to " + getParserName() + " should be a single value if element-wise is being used.");
        }

        ProductStatistic productStatistic = new ProductStatistic(name, elementwise, constants);

        for (int i = 0; i < xo.getChildCount(); i++) {
            Object child = xo.getChild(i);
            if (child instanceof Statistic) {
                try {
                    Statistic statistic = (Statistic)child;

                    if (constants != null && constants.length != 1 && constants.length != statistic.getDimension()) {
                        throw new XMLParseException("The constants given to " + getParserName() + " is not of the same dimension as the statistics");
                    }

                    productStatistic.addStatistic(statistic);
                } catch (IllegalArgumentException iae) {
                    throw new XMLParseException("Statistic added to " + getParserName() + " element is not of the same dimension");
                }
            } else {
                throw new XMLParseException("Unknown element found in " + getParserName() + " element:" + child);
            }
        }

        return productStatistic;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element returns a statistic that is the product of the child statistics.";
    }

    public Class getReturnType() {
        return ProductStatistic.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newBooleanRule(ELEMENT_WISE, true),
            AttributeRule.newStringRule(Statistic.NAME, true),
            AttributeRule.newDoubleArrayRule(CONSTANT, true, "A constant to be multiplied with the parameters. If element-wise then this can be a vector of constants."),
            new ElementRule(Statistic.class, 1, Integer.MAX_VALUE)
    };
}

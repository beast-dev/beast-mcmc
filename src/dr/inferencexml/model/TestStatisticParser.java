/*
 * TestStatisticParser.java
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

package dr.inferencexml.model;

import dr.inference.model.Statistic;
import dr.inference.model.TestStatistic;
import dr.util.Attribute;
import dr.xml.*;

/**
 */
public class TestStatisticParser extends AbstractXMLObjectParser {

    public static String TEST_STATISTIC = "testStatistic";
    public static String TEST = "test";

    private static final String SEQUALS = "equals";
    private static final String SGREATER_THAN = "greaterThan";
    private static final String SLESS_THAN = "lessThan";
    private static final String SINSIDE = "inside";
    private static final String SOUTSIDE = "outside";

    public String getParserName() {
        return TEST_STATISTIC;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        String name = xo.getAttribute(Statistic.NAME, xo.hasId() ? xo.getId() : "");
        Attribute attr = (Attribute) xo.getChild(Attribute.class);
        double testValue1;

        TestStatistic statistic;

        if (xo.hasChildNamed(SEQUALS)) {
            Attribute attr2 = (Attribute) xo.getElementFirstChild(SEQUALS);
            statistic = new TestStatistic(name, attr, attr2, TestStatistic.EQUALS);
        } else if (xo.hasChildNamed(SGREATER_THAN)) {
            Attribute attr2 = (Attribute) xo.getElementFirstChild(SGREATER_THAN);
            statistic = new TestStatistic(name, attr, attr2, TestStatistic.GREATER_THAN);
        } else if (xo.hasChildNamed(SLESS_THAN)) {
            Attribute attr2 = (Attribute) xo.getElementFirstChild(SLESS_THAN);
            statistic = new TestStatistic(name, attr, attr2, TestStatistic.LESS_THAN);
        } else if (xo.hasAttribute(SEQUALS)) {
            testValue1 = xo.getDoubleAttribute(SEQUALS);
            statistic = new TestStatistic(name, attr, testValue1, TestStatistic.EQUALS);
        } else if (xo.hasAttribute(SGREATER_THAN)) {
            testValue1 = xo.getDoubleAttribute(SGREATER_THAN);
            statistic = new TestStatistic(name, attr, testValue1, TestStatistic.GREATER_THAN);
        } else if (xo.hasAttribute(SLESS_THAN)) {
            testValue1 = xo.getDoubleAttribute(SLESS_THAN);
            statistic = new TestStatistic(name, attr, testValue1, TestStatistic.LESS_THAN);
        } else if (xo.hasAttribute(SINSIDE)) {
            double[] values = xo.getDoubleArrayAttribute(SINSIDE);
            if (values.length != 2)
                throw new XMLParseException("inside attribute of test element requires two values");
            statistic = new TestStatistic(name, attr, values[0], values[1], TestStatistic.INSIDE);
        } else if (xo.hasAttribute(SOUTSIDE)) {
            double[] values = xo.getDoubleArrayAttribute(SOUTSIDE);
            if (values.length != 2)
                throw new XMLParseException("outside attribute of test element requires two values");
            statistic = new TestStatistic(name, attr, values[0], values[1], TestStatistic.OUTSIDE);
        } else throw new XMLParseException();

        return statistic;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element represents a boolean statistic that returns 1 " +
                        "if the conditions are met and 0 otherwise.";
    }

    public Class getReturnType() {
        return TestStatistic.class;
    }

    @Override
    public String[] getParserNames() {
        return new String[]{getParserName(), TEST};
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new StringAttributeRule("name", "A name for this statistic, for logging purposes", true),
            new ElementRule(Attribute.class),
            new XORRule(
                    new XMLSyntaxRule[]{
                            new ElementRule(SEQUALS, Attribute.class),
                            new ElementRule(SGREATER_THAN, Attribute.class),
                            new ElementRule(SLESS_THAN, Attribute.class),
                            AttributeRule.newDoubleRule(SEQUALS),
                            AttributeRule.newDoubleRule(SGREATER_THAN),
                            AttributeRule.newDoubleRule(SLESS_THAN),
                            AttributeRule.newDoubleArrayRule(SINSIDE),
                            AttributeRule.newDoubleArrayRule(SOUTSIDE)
                    }
            )
    };
}

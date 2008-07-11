/*
 * TestStatistic.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
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

package dr.inference.model;

import dr.util.Attribute;
import dr.xml.*;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @version $Id: TestStatistic.java,v 1.12 2005/07/11 14:06:25 rambaut Exp $
 */
public class TestStatistic extends BooleanStatistic {

    public static String TEST_STATISTIC = "test";
    public static String NAME = "name";

    private Attribute attribute = null;
    private Attribute attribute2 = null;
    private double testValue1, testValue2;
    private int mode;

    private static final int EQUALS = 0;
    private static final int GREATER_THAN = 1;
    private static final int LESS_THAN = 2;
    private static final int INSIDE = 3;
    private static final int OUTSIDE = 4;


    public TestStatistic(String name, Attribute attr, double value, int mode) {
        super(name);

        attribute = attr;
        testValue1 = value;
        this.mode = mode;
    }

    /**
     * A constructor for modes that require 2 values (i.e., INSIDE and OUTSIDE).
     */
    public TestStatistic(String name, Attribute attr, double value1, double value2, int mode) {
        super(name);

        attribute = attr;
        testValue1 = value1;
        testValue2 = value2;
        this.mode = mode;
    }

    /**
     * A constructor that compares 2 attributes (i.e., GREATER_THAN or LESS_THAN).
     */
    public TestStatistic(String name, Attribute attr1, Attribute attr2, int mode) {
        super(name);

        attribute = attr1;
        attribute2 = attr2;
        this.mode = mode;
    }

    public int getDimension() {
        return 1;
    }

    /**
     * @return boolean result of test.
     */
    public boolean getBoolean(int i) {

        double num;
        if (attribute instanceof Statistic) {
            num = ((Statistic) attribute).getStatisticValue(0);
        } else {
            num = ((Number) attribute.getAttributeValue()).doubleValue();
        }

        if (attribute2 != null) {
            if (attribute2 instanceof Statistic) {
                testValue1 = ((Statistic) attribute2).getStatisticValue(0);
            } else {
                testValue1 = ((Number) attribute2.getAttributeValue()).doubleValue();
            }
        }

        switch (mode) {
            case EQUALS:
                if (num == testValue1) return true;
                break;
            case GREATER_THAN:
                if (num > testValue1) return true;
                break;
            case LESS_THAN:
                if (num < testValue1) return true;
                break;
            case INSIDE:
                if (num > testValue1 && num < testValue2) return true;
                break;
            case OUTSIDE:
                if (num < testValue1 && num > testValue2) return true;
                break;
        }
        return false;
    }


    public static dr.xml.XMLObjectParser PARSER = new dr.xml.AbstractXMLObjectParser() {

        public String getParserName() {
            return TEST_STATISTIC;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            String name = xo.getAttribute(NAME, xo.hasId() ? xo.getId() : "");
            Attribute attr = (Attribute) xo.getChild(Attribute.class);
            double testValue1;

            TestStatistic statistic;

            if (xo.hasChildNamed("equals")) {
                Attribute attr2 = (Attribute) xo.getElementFirstChild("equals");
                statistic = new TestStatistic(name, attr, attr2, EQUALS);
            } else if (xo.hasChildNamed("greaterThan")) {
                Attribute attr2 = (Attribute) xo.getElementFirstChild("greaterThan");
                statistic = new TestStatistic(name, attr, attr2, GREATER_THAN);
            } else if (xo.hasChildNamed("lessThan")) {
                Attribute attr2 = (Attribute) xo.getElementFirstChild("lessThan");
                statistic = new TestStatistic(name, attr, attr2, LESS_THAN);
            } else if (xo.hasAttribute("equals")) {
                testValue1 = xo.getDoubleAttribute("equals");
                statistic = new TestStatistic(name, attr, testValue1, EQUALS);
            } else if (xo.hasAttribute("greaterThan")) {
                testValue1 = xo.getDoubleAttribute("greaterThan");
                statistic = new TestStatistic(name, attr, testValue1, GREATER_THAN);
            } else if (xo.hasAttribute("lessThan")) {
                testValue1 = xo.getDoubleAttribute("lessThan");
                statistic = new TestStatistic(name, attr, testValue1, LESS_THAN);
            } else if (xo.hasAttribute("inside")) {
                double[] values = xo.getDoubleArrayAttribute("inside");
                if (values.length != 2)
                    throw new XMLParseException("inside attribute of test element requires two values");
                statistic = new TestStatistic(name, attr, values[0], values[1], INSIDE);
            } else if (xo.hasAttribute("outside")) {
                double[] values = xo.getDoubleArrayAttribute("outside");
                if (values.length != 2)
                    throw new XMLParseException("outside attribute of test element requires two values");
                statistic = new TestStatistic(name, attr, values[0], values[1], OUTSIDE);
            } else throw new XMLParseException();

            return statistic;
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return
                    "This element represents a boolean statistic that returns 1 " +
                            "if the conditions are met and 0 otherwise.";
        }

        public Class getReturnType() {
            return BooleanStatistic.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new StringAttributeRule("name", "A name for this statistic, for logging purposes", true),
                new ElementRule(Attribute.class),
                new XORRule(
                        new XMLSyntaxRule[]{
                                new ElementRule("equals", Attribute.class),
                                new ElementRule("greaterThan", Attribute.class),
                                new ElementRule("lessThan", Attribute.class),
                                AttributeRule.newDoubleRule("equals"),
                                AttributeRule.newDoubleRule("greaterThan"),
                                AttributeRule.newDoubleRule("lessThan"),
                                AttributeRule.newDoubleArrayRule("inside"),
                                AttributeRule.newDoubleArrayRule("outside")
                        }
                )
		};
	
	};
}

/*
 * BeastUnitTest.java
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

package dr.xml.unittest;

import dr.xml.*;

import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * @author Gabriel Hassler
 * @author Marc A. Suchard
 * @author Xinghua Tao
 */


public class BeastUnitTest implements Reportable {

    private final String message;
    private final String[] actual;
    private final String expected;
    private final int[] indices;
    private Boolean pass;
    private final AssertType assertType;
    private final boolean equalMode;

    public BeastUnitTest(String message, String[] actual, String expected,
                         AssertType assertType, int[] indices, boolean equalMode) {
        this.message = message;
        this.actual = actual;
        this.expected = expected;
        this.assertType = assertType;
        this.indices = indices;
        this.equalMode = equalMode;
    }

    public void execute() {
        for (int i = 0; i < actual.length; i++) {
            boolean result = equalMode
                    ? assertType.equivalent(actual[i], indices, expected)
                    : !assertType.equivalent(actual[i], indices, expected);
            if (!result) {
                failCheck(i);
            }
        }
        pass = true;
    }

    private void failCheck(int i) {
        String string = formatName()
                + ": '" + actual[i] + "' != '" + expected + "'";
        System.err.println(string);
        System.exit(-1);
    }

    private String formatName() {
        return "assert" + ((message != null) ? (" " + message) : "");
    }

    private String getPass() {
        return pass != null ?
                (pass ? "passed" : "fail") :
                "not executed";
    }

    @Override
    public String getReport() {
        return formatName() + ": " + getPass();
    }

    interface AssertType {

        boolean equivalent(String a, int[] aIndices, String b);

        class StringAssert implements AssertType {
            @Override
            public boolean equivalent(String a, int[] aIndices, String b) {
                return a.compareTo(b) == 0;
            }
        }

        class DoubleAssert implements AssertType {

            private final double tolerance;
            private final String stripChars;
            private final ToleranceType toleranceType;

            DoubleAssert(double tolerance, ToleranceType tolType, String stripChars) {
                this.tolerance = tolerance;
                this.stripChars = stripChars;
                this.toleranceType = tolType;
            }

            @Override
            public boolean equivalent(String a, int[] aIndices, String b) {

                double[] lhs = parseArray(a, aIndices);
                double[] rhs = parseArray(b, null);
                if (lhs.length != rhs.length) {
                    System.err.println("The dimensions of the \"actual\" and \"expected\" values are not the same.");
                    return false;
                }

                for (int i = 0; i < lhs.length; ++i) {
                    if (!toleranceType.close(lhs[i], rhs[i], tolerance)) {
                        System.err.println("Dimension " + (i + 1) + " of \"actual\" does not match \" expected\". (" + lhs[i] + " != " + rhs[i] + ")");
                        return false;
                    }
                }

                return true;
            }

            private double[] parseArray(String string, int[] indices) {
                string = string.replaceAll(",", " ");
                string = string.replaceAll("[" + stripChars + "]", " ");
                string = string.trim();
                String[] strings = string.split("\\s+");
                final double[] reals;
                if (indices == null) {
                    reals = new double[strings.length];

                    for (int i = 0; i < strings.length; ++i) {
                        reals[i] = Double.parseDouble(strings[i]);
                    }
                } else {
                    reals = new double[indices.length];

                    int dim = 0;
                    for (int i : indices) {
                        reals[dim] = Double.parseDouble(strings[i]);
                        dim++;
                    }
                }

                return reals;
            }

            enum ToleranceType {
                RELATIVE {
                    @Override
                    boolean close(double lhs, double rhs, double tolerance) {
                        double tol = Math.abs(tolerance * rhs);
                        return ToleranceType.ABSOLUTE.close(lhs, rhs, tol);
                    }
                },
                ABSOLUTE {
                    @Override
                    boolean close(double lhs, double rhs, double tolerance) {
                        return Math.abs(lhs - rhs) < tolerance;
                    }
                };

                abstract boolean close(double lhs, double rhs, double tolerance);
            }
        }
    }

    private static final String CHECK = "assertEqual";
    private static final String EQUAL = "equal";
    private static final String MESSAGE = "message";
    private static final String EXPECTED = "expected";
    private static final String ACTUAL = "actual";
    private static final String REGEX = "regex";
    private static final String TOLERANCE_STRING = "tolerance";
    private static final String VERBOSE = "verbose";
    private static final String STRIP_CHARACTERS = "charactersToStrip";
    private static final String TOLERANCE_TYPE = "toleranceType";
    private static final String ABSOLUTE = "absolute";
    private static final String RELATIVE = "relative";
    private static final String INDICES = "actualIndices";

    public static AbstractXMLObjectParser PARSER = new AbstractXMLObjectParser() {

        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            String message = null;
            if (xo.hasChildNamed(MESSAGE)) {
                message = xo.getChild(MESSAGE).getStringChild(0);
            }

            boolean equal = xo.getAttribute(EQUAL, true);

            String[] expectedArray = parseValues(xo.getChild(EXPECTED));
            if (expectedArray.length != 1) {
                throw new XMLParseException("There should only be one '" + EXPECTED + "' value.");
            }
            String expected = expectedArray[0];

            String[] actual = parseValues(xo.getChild(ACTUAL));
            String stripChars = xo.getAttribute(STRIP_CHARACTERS, ",");

            AssertType assertType;

            if (xo.hasAttribute(TOLERANCE_STRING)) {
                double tolerance = xo.getDoubleAttribute(TOLERANCE_STRING);
                String tolTypeString = xo.getAttribute(TOLERANCE_TYPE, ABSOLUTE);

                AssertType.DoubleAssert.ToleranceType tolType;

                if (tolTypeString.equalsIgnoreCase(ABSOLUTE)) {
                    tolType = AssertType.DoubleAssert.ToleranceType.ABSOLUTE;

                } else if (tolTypeString.equalsIgnoreCase(RELATIVE)) {
                    tolType = AssertType.DoubleAssert.ToleranceType.RELATIVE;

                } else {
                    throw new XMLParseException("The optional attribute " + TOLERANCE_TYPE + " must be either \"" +
                            RELATIVE + "\" or \"" + ABSOLUTE + "\"");
                }

                assertType = new AssertType.DoubleAssert(tolerance, tolType, stripChars);

            } else {
                assertType = new AssertType.StringAssert();
            }
            int[] indices = null;
            if (xo.hasAttribute(INDICES)) indices = xo.getIntegerArrayAttribute(INDICES);

            BeastUnitTest unitTest = new BeastUnitTest(message, actual, expected, assertType, indices, equal);
            unitTest.execute();

            if (xo.getAttribute(VERBOSE, false)) {
                Logger.getLogger("dr.xml.unittest").info(unitTest.getReport());
            }

            return unitTest;
        }

        private String[] parseValues(XMLObject xo) throws XMLParseException {
            int nChildren = xo.getChildCount();
            String[] rawStrings = new String[nChildren];

            for (int i = 0; i < nChildren; i++) {
                if (xo.getChild(i) instanceof Reportable) {
                    Reportable reportable = (Reportable) xo.getChild(i);
                    rawStrings[i] = reportable.getReport();
                } else {
                    rawStrings[i] = xo.getStringChild(i);
                }
            }

            if (xo.hasAttribute(REGEX)) {
                Pattern pattern = Pattern.compile(xo.getStringAttribute(REGEX));
                for (int i = 0; i < nChildren; i++) {
                    Matcher matcher = pattern.matcher(rawStrings[i]);
                    if (matcher.find()) {
                        rawStrings[i] = matcher.group(1);
                    }
                }
            }

            return rawStrings;
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        @Override
        public String getParserDescription() {
            return null;
        }

        @Override
        public Class getReturnType() {
            return BeastUnitTest.class;
        }

        @Override
        public String getParserName() {
            return CHECK;
        }
    };

    private final static XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            AttributeRule.newBooleanRule(EQUAL, true),
            new ElementRule(EXPECTED, new XMLSyntaxRule[]{
                    new XORRule(
                            new ElementRule(Reportable.class),
                            new ElementRule(String.class)
                    ),
                    AttributeRule.newStringRule(REGEX, true),
            }),
            new ElementRule(ACTUAL, new XMLSyntaxRule[]{
                    new XORRule(
                            new ElementRule(Reportable.class, 1, Integer.MAX_VALUE),
                            new ElementRule(String.class, 1, Integer.MAX_VALUE)
                    ),
                    AttributeRule.newStringRule(REGEX, true),
            }),
            new ElementRule(MESSAGE, new XMLSyntaxRule[]{
                    new ElementRule(String.class),
            }, true),
            AttributeRule.newDoubleRule(TOLERANCE_STRING, true),
            AttributeRule.newBooleanRule(VERBOSE, true),
            AttributeRule.newStringRule(TOLERANCE_TYPE, true),
            AttributeRule.newIntegerArrayRule(INDICES, true)
    };
}

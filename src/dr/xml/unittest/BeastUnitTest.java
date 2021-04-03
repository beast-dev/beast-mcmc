package dr.xml.unittest;

import dr.xml.*;

import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * @author Gabriel Hassler
 * @author Marc A. Suchard
 */

public class BeastUnitTest implements Reportable {

    private final String message;
    private final String[] actual;
    private final String expected;

    private Boolean pass;

    private final AssertType assertType;

    public BeastUnitTest(String message, String[] actual, String expected,
                         AssertType assertType) {
        this.message = message;
        this.actual = actual;
        this.expected = expected;
        this.assertType = assertType;
    }

    public void execute() {
        for (int i = 0; i < actual.length; i++) {
            if (!assertType.equivalent(actual[i], expected)) {
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

        boolean equivalent(String a, String b);

        class StringAssert implements AssertType {

            @Override
            public boolean equivalent(String a, String b) {
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
            public boolean equivalent(String a, String b) {

                double[] lhs = parseArray(a);
                double[] rhs = parseArray(b);

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

            private double[] parseArray(String string) {
                string = string.replaceAll(",", " ");
                string = string.replaceAll("[" + stripChars + "]", " ");
                string = string.trim();
                String[] strings = string.split("\\s+");
                double[] reals = new double[strings.length];
                for (int i = 0; i < strings.length; ++i) {
                    reals[i] = Double.valueOf(strings[i]);
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

    public static AbstractXMLObjectParser PARSER = new AbstractXMLObjectParser() {

        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            String message = null;
            if (xo.hasChildNamed(MESSAGE)) {
                message = xo.getChild(MESSAGE).getStringChild(0);
            }

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

            BeastUnitTest unitTest = new BeastUnitTest(message, actual, expected, assertType);
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
            AttributeRule.newStringRule(TOLERANCE_TYPE, true)
    };
}

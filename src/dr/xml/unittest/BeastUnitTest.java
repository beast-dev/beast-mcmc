package dr.xml.unittest;

import dr.xml.*;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * @author Gabriel Hassler
 * @author Marc A. Suchard
 */

public class BeastUnitTest {

    private static final double TOLERANCE = 1e-6;

    private final String message;
    private final String actual;
    private final String expected;
    private final double tolerance;

    private final AssertType assertType;

    public BeastUnitTest(String message, String actual, String expected) {
        this(message, actual, expected, TOLERANCE);
    }

    public BeastUnitTest(String message, String actual, String expected, double tolerance) {
        this.message = message;
        this.actual = actual;
        this.expected = expected;
        this.tolerance = tolerance;

        assertType = AssertType.STRING;
    }

    public void execute() {
        if (!assertType.equivalent(actual, expected, tolerance)) {
            failCheck();
        }
    }

    private void failCheck() {
        String string = "assert " + ((message != null) ? message : "")
                + actual + " != " + expected;
        System.err.println(string);
        System.exit(-1);
    }

    enum AssertType {
        STRING {
            @Override
            boolean equivalent(String a, String b, double tolerance) {
                return a.compareTo(b) == 0;
            }
        },
        DOUBLE {
            @Override
            boolean equivalent(String a, String b, double tolerance) {
                double lhs = Double.valueOf(a);
                double rhs = Double.valueOf(b);
                return close(lhs, rhs, tolerance);
            }
        },
        DOUBLE_ARRAY {
            @Override
            boolean equivalent(String a, String b, double tolerance) {

                double[] lhs = parseArray(a);
                double[] rhs = parseArray(b);

                if (lhs.length != rhs.length) {
                    return false;
                }

                for (int i = 0; i < lhs.length; ++i) {
                    if (!close(lhs[i], rhs[i], tolerance)) {
                        return false;
                    }
                }

                return true;
            }

            double[] parseArray(String string) {
                string = string.replaceAll(",", " ");
                String[] strings = string.split("\\s+");
                double[] reals = new double[string.length()];
                for (int i = 0; i < strings.length; ++i) {
                    reals[i] = Double.valueOf(strings[i]);
                }

                return reals;
            }
        };

        abstract boolean equivalent(String a, String b, double tolerance);

        private static boolean close(double lhs, double rhs, double tolerance) {
            return Math.abs(lhs - rhs) < tolerance;
        }
    }

    private static final String CHECK = "assertEqual";
    private static final String MESSAGE = "message";
    private static final String EXPECTED = "expected";
    private static final String ACTUAL = "actual";
    private static final String REGEX = "regex";

    public static AbstractXMLObjectParser PARSER = new AbstractXMLObjectParser() {

        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            String message = null;
            if (xo.hasChildNamed(MESSAGE)) {
                message = xo.getChild(MESSAGE).getStringChild(0);
            }

            String expected = parseValue(xo.getChild(EXPECTED));
            String actual = parseValue(xo.getChild(ACTUAL));

            BeastUnitTest unitTest = new BeastUnitTest(message, actual, expected);
            unitTest.execute();
            
            return unitTest;
        }

        private String parseValue(XMLObject xo) throws XMLParseException {

            String rawString;

            if (xo.getChild(0) instanceof Reportable) {
                Reportable reportable = (Reportable) xo.getChild(0);
                rawString = reportable.getReport();
            } else {
                rawString = xo.getStringChild(0);
            }

            if (xo.hasAttribute(REGEX)) {
                Pattern pattern = Pattern.compile(xo.getStringAttribute(REGEX));
                Matcher matcher = pattern.matcher(rawString);
                if (matcher.find()) {
                    rawString = matcher.group(1);
                }
            }

            return rawString;
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
                            new ElementRule(Reportable.class),
                            new ElementRule(String.class)
                    ),
                    AttributeRule.newStringRule(REGEX, true),
            }),
            new ElementRule(MESSAGE, new XMLSyntaxRule[]{
                    new ElementRule(String.class),

            }, true),
    };
}



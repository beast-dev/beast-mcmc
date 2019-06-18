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

    private static final double TOLERANCE = 1e-6;

    private final String message;
    private final String actual;
    private final String expected;
    private final double tolerance;

    private Boolean pass;

    private final AssertType assertType;

    public BeastUnitTest(String message, String actual, String expected) {
        this(message, actual, expected,
                AssertType.STRING, TOLERANCE);
    }

    public BeastUnitTest(String message, String actual, String expected,
                         AssertType assertType,
                         double tolerance) {
        this.message = message;
        this.actual = actual;
        this.expected = expected;
        this.assertType = assertType;
        this.tolerance = tolerance;
    }

    public void execute() {
        if (!assertType.equivalent(actual, expected, tolerance)) {
            failCheck();
        }

        pass = true;
    }

    private void failCheck() {
        String string = formatName()
                + ": '" + actual + "' != '" + expected + "'";
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

            private double[] parseArray(String string) {
                string = string.replaceAll(",", " ");
                String[] strings = string.split("\\s+");
                double[] reals = new double[strings.length];
                for (int i = 0; i < strings.length; ++i) {
                    reals[i] = Double.valueOf(strings[i]);
                }

                return reals;
            }

            private boolean close(double lhs, double rhs, double tolerance) {
                return Math.abs(lhs - rhs) < tolerance;
            }
        };

        abstract boolean equivalent(String a, String b, double tolerance);
    }

    private static final String CHECK = "assertEqual";
    private static final String MESSAGE = "message";
    private static final String EXPECTED = "expected";
    private static final String ACTUAL = "actual";
    private static final String REGEX = "regex";
    private static final String TOLERANCE_STRING = "tolerance";
    private static final String VERBOSE = "verbose";

    public static AbstractXMLObjectParser PARSER = new AbstractXMLObjectParser() {

        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            String message = null;
            if (xo.hasChildNamed(MESSAGE)) {
                message = xo.getChild(MESSAGE).getStringChild(0);
            }

            String expected = parseValue(xo.getChild(EXPECTED));
            String actual = parseValue(xo.getChild(ACTUAL));

            BeastUnitTest unitTest;
            if (xo.hasAttribute(TOLERANCE_STRING)) {
                double tolerance = xo.getAttribute(TOLERANCE_STRING, TOLERANCE);
                unitTest = new BeastUnitTest(message, actual, expected,
                        AssertType.DOUBLE, tolerance);
            } else {
                unitTest = new BeastUnitTest(message, actual, expected);
            }
            unitTest.execute();

            if (xo.getAttribute(VERBOSE, false)) {
                Logger.getLogger("dr.xml.unittest").info(unitTest.getReport());
            }
            
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
            AttributeRule.newDoubleRule(TOLERANCE_STRING, true),
            AttributeRule.newBooleanRule(VERBOSE, true),
    };
}



package dr.xml;

import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class Check {

    private static final Pattern doublePattern = Pattern.compile("[-+]?\\d+(\\.\\d+)?");
    private final Report report;
    private static final double tolerance = 1e-6;


    public Check(Report report, String[] regexes, String[] values) {

        this.report = report;
        String stringReport = getReportAsString();

        for (int i = 0; i < regexes.length; i++) {
            String stringValue = parseDoubleRegex(stringReport, regexes[i]);
//            boolean passes = stringValue.equals(values[i]);
            boolean passes = checkEqualsDouble(stringValue, values[i]);
            if (!passes) {
                failCheck(regexes[i], stringValue, values[i]);
            }
        }

    }

    private void failCheck(String name, String reportValue, String trueValue) {
        System.out.println("Report returned " + reportValue + " for " + name + "."
                + " The true value is " + trueValue + ".");
        System.exit(-1);
    }

    private boolean checkEqualsDouble(String s1, String s2) {
        double d1 = Double.parseDouble(s1);
        double d2 = Double.parseDouble(s2);
        return Math.abs(d1 - d2) > tolerance ? false : true;

    }

    private String getReportAsString() {

        StringWriter sWriter = new StringWriter();
        PrintWriter pWriter = new PrintWriter(sWriter);
        report.setOutput(pWriter);

        report.createReport();

        return sWriter.toString();

    }

    private String parseDoubleRegex(String report, String matchBase) {

        Pattern keyPattern = stringMatchPattern(matchBase);
        Matcher keyMatcher = keyPattern.matcher(report);
        keyMatcher.find();

        if (keyMatcher.end() == -1) {
            throw new RuntimeException("Did not find \"" + matchBase + "\" in the report.");
        }

        Matcher matcher = doublePattern.matcher(report);
        matcher.find(keyMatcher.end());

        int sInd = matcher.start();
        int lInd = matcher.end();
        int stringLength = lInd - sInd;

        char[] doubleChars = new char[stringLength];
        report.getChars(sInd, lInd, doubleChars, 0);
        String doubleString = new String(doubleChars);

        return doubleString;
    }

    private Pattern stringMatchPattern(String baseMatch) {
        String regex = "\\s+" + baseMatch + "[\\s+|:]";
        return Pattern.compile(regex);
    }


    private static final String CHECK = "check";
    private static final String REGEX = "regex";
    private static final String MATCH = "match";


    public static AbstractXMLObjectParser PARSER = new AbstractXMLObjectParser() {
        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            Report report = (Report) xo.getChild(Report.class);


            List<XMLObject> matches = xo.getAllChildren(MATCH);
            int nMatch = matches.size();
            String[] regexes = new String[nMatch];
            String[] values = new String[nMatch];
            for (int i = 0; i < nMatch; i++) {
                XMLObject match = matches.get(i);
                regexes[i] = (String) match.getAttribute(REGEX);
                values[i] = match.getStringChild(0);
            }

            return new Check(report, regexes, values);
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
            return null;
        }

        @Override
        public String getParserName() {
            return CHECK;
        }
    };

    private final static XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(Report.class),
            new ElementRule(MATCH, new XMLSyntaxRule[]{
                    AttributeRule.newStringRule(REGEX)
            }, 1, Integer.MAX_VALUE)};
};



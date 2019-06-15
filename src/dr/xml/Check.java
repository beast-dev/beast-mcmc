package dr.xml;

import java.io.*;
import java.util.ArrayList;
import java.util.List;

public class Check {



    public Check(Report report, String[] regexes, String[] values) {
        String stringReport = getReportAsString(report);
        //TODO: parse regexes & values
        //TODO: make comparisons
    }

    private String getReportAsString(Report report){

        StringWriter sWriter = new StringWriter();
        PrintWriter pWriter = new PrintWriter(sWriter);
        report.setOutput(pWriter);

        report.createReport();

        return sWriter.toString();

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
            for (int i = 0; i < nMatch; i++){
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



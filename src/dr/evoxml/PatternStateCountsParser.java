package dr.evoxml;

import dr.evolution.alignment.PatternList;
import dr.evolution.alignment.PatternStateCounts;
import dr.util.FileHelpers;
import dr.xml.*;

import java.io.PrintWriter;

import static dr.inferencexml.loggers.LoggerParser.getLogFile;

public class PatternStateCountsParser extends AbstractXMLObjectParser {

    public static final String PATTERN_STATE_COUNTS = "patternStateCounts";
    private static final String FILE_NAME = FileHelpers.FILE_NAME;

    @Override
    public String getParserName() {
        return PATTERN_STATE_COUNTS;
    }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        PatternList patternList = (PatternList) xo.getChild(PatternList.class);
        PrintWriter pw = getLogFile(xo, getParserName());

        new PatternStateCounts(patternList).report(pw);

        pw.flush();
        if (xo.hasAttribute(FILE_NAME)) {
            pw.close();
        }
        return null;
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newStringRule(FILE_NAME, true),
            new ElementRule(PatternList.class),
    };

    @Override
    public String getParserDescription() {
        return "For each unique site pattern, prints counts, proportions, entropy, distinct-state count, and gap/ambiguous count.";
    }

    @Override
    public Class getReturnType() {
        return Object.class;
    }
}

package dr.evomodelxml.tree;

import java.io.PrintWriter;

import dr.evomodel.tree.HiddenLinkageLogger;
import dr.evomodel.tree.HiddenLinkageModel;
import dr.inference.loggers.LogFormatter;
import dr.inference.loggers.TabDelimitedFormatter;
import dr.inferencexml.loggers.LoggerParser;
import dr.xml.AttributeRule;
import dr.xml.ElementRule;
import dr.xml.StringAttributeRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

/**
 * @author Aaron Darling
 */
public class HiddenLinkageLoggerParser extends LoggerParser {
    public static final String LOG_HIDDEN_LINKAGE = "logHiddenLinkage";

    public String getParserName() {
        return LOG_HIDDEN_LINKAGE;
    }

    /**
     * @return an object based on the XML element it was passed.
     */
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
    	HiddenLinkageModel hlm = (HiddenLinkageModel)xo.getChild(HiddenLinkageModel.class);
        // logEvery of zero only displays at the end
        final int logEvery = xo.getAttribute(LOG_EVERY, 0);
        final PrintWriter pw = getLogFile(xo, getParserName());
        final LogFormatter formatter = new TabDelimitedFormatter(pw);

        return new HiddenLinkageLogger(hlm, formatter, logEvery);
    }

    public String getParserDescription() {
        return "Logs a linkage groups for metagenomic reads to a file";
    }

    public Class getReturnType() {
        return HiddenLinkageLogger.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }
    private final XMLSyntaxRule[] rules = {
            AttributeRule.newIntegerRule(LOG_EVERY, true),
            new StringAttributeRule(FILE_NAME,
                    "The name of the file to send log output to. " +
                            "If no file name is specified then log is sent to standard output", true),
            new ElementRule(HiddenLinkageModel.class, "The linkage model which is to be logged"),
    };
}

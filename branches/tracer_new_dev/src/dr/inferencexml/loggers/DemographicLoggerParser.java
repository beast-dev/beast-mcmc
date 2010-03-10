package dr.inferencexml.loggers;

import dr.evomodel.coalescent.DemographicLogger;
import dr.evomodel.coalescent.DemographicReconstructor;
import dr.inference.loggers.LogFormatter;
import dr.inference.loggers.TabDelimitedFormatter;
import dr.xml.*;

import java.io.PrintWriter;

/**
 * Parses an element from an DOM document into a ExponentialGrowth.
 */
public class DemographicLoggerParser extends AbstractXMLObjectParser {
    public static final String DEMOGRAPHIC_LOG = "demographicLog";
    public static final String TITLE = "title";
    public static final String FILE_NAME = "fileName";
    public static final String LOG_EVERY = "logEvery";

    public String getParserName() {
        return DEMOGRAPHIC_LOG;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        final int logEvery = xo.getIntegerAttribute(LOG_EVERY);

        final PrintWriter pw = LoggerParser.getLogFile(xo, getParserName());

        final LogFormatter formatter = new TabDelimitedFormatter(pw);

        DemographicReconstructor reconstructor = (DemographicReconstructor) xo.getChild(DemographicReconstructor.class);

        return new DemographicLogger(reconstructor, formatter, logEvery);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A demographic model of constant population size followed by logistic growth.";
    }

    public Class getReturnType() {
        return DemographicLogger.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            AttributeRule.newIntegerRule(LOG_EVERY),
            new StringAttributeRule(FILE_NAME,
                    "The name of the file to send log output to. " +
                            "If no file name is specified then log is sent to standard output", true),
            new ElementRule(DemographicReconstructor.class)
    };
}

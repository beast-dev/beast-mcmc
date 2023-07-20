package dr.inferencexml.loggers;

import dr.inference.loggers.TimeLogger;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

public class TimeLoggerXMLParser extends AbstractXMLObjectParser {
    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        return new TimeLogger();
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[0];
    }

    @Override
    public String getParserDescription() {
        return "saves elapsed seconds to log file";
    }

    @Override
    public Class getReturnType() {
        return TimeLogger.class;
    }

    @Override
    public String getParserName() {
        return "timeLogger";
    }
}

package dr.inferencexml.loggers;

import dr.inference.loggers.PrngStateLogger;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

public class PrngStateLoggerParser extends AbstractXMLObjectParser {

    private static final String PARSER_NAME = "prngStateLogger";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        return new PrngStateLogger();
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[0];
    }

    @Override
    public String getParserDescription() {
        return "Logs the internal PRNG state";
    }

    @Override
    public Class getReturnType() {
        return PrngStateLogger.class;
    }

    @Override
    public String getParserName() { return PARSER_NAME; }
}

package @PLUGIN_FULL_PACKAGE@;

@PLUGIN_IMPORTS@
import dr.xml.*;

import java.util.logging.Logger;

public class @PLUGIN_CLASS@Parser extends AbstractXMLObjectParser {

    public String getParserName() {
    	return "@PLUGIN_XML_ELEMENT@";
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
    	@PLUGIN_XML_PARSER_STUB@
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
    	return "This element represents an instance that comes from the @PLUGIN_CLASS@";
    }

    public Class getReturnType() {
    	return @PLUGIN_PARSER_RETURN_TYPE@;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
    	return rules;
    }

    private final XMLSyntaxRule[] rules = {
    	@PLUGIN_XML_SYNTAX_RULES@
    };
}

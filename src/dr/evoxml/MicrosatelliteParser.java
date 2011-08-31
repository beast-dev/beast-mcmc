package dr.evoxml;

import dr.xml.*;
import dr.evolution.datatype.Microsatellite;


/**
 * @author Chieh-Hsi Wu
 *
 * Microsatellite data type parser
 *
 */
public class MicrosatelliteParser extends AbstractXMLObjectParser {

    public static final String MICROSAT = "microsatellite";

    public String getParserName() {
        return MICROSAT;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        int min = xo.getIntegerAttribute("min");
        int max = xo.getIntegerAttribute("max");
        int unitLength = xo.getIntegerAttribute("unitLength");

        return new Microsatellite(min, max, unitLength);
    }

    public String getParserDescription() {
        return "This element represents a microsatellite data type.";
    }

    public String getExample() {
        return "<microsatellite min=\"0\" max=\"20\" unitLength=\"2\"/>";
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                new AndRule(new XMLSyntaxRule[]{
                    AttributeRule.newIntegerRule("min"),
                    AttributeRule.newIntegerRule("max"),
                    AttributeRule.newIntegerRule("unitLength")})
        };
    }


    public Class getReturnType() {
        return Microsatellite.class;
    }
}
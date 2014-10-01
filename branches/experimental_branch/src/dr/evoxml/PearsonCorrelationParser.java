package dr.evoxml;

import dr.xml.*;
import dr.inference.model.Parameter;
import dr.inference.model.PearsonCorrelation;

/**
 * @author Simon Greenhill
 */

public class PearsonCorrelationParser extends AbstractXMLObjectParser {

    public static final String PEARSON_CORRELATION = "pearsonCorrelation";
    public static final String LOG = "log";

    public String getParserName() { return PEARSON_CORRELATION; }
                                             
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        boolean log = xo.getAttribute(LOG, false);
        Parameter X = (Parameter)xo.getChild(0);
        Parameter Y = (Parameter)xo.getChild(1);

        // System.out.println("Correlating " + X + " with " + Y + " using log = " + log);
        PearsonCorrelation pearsonCorrelation = new PearsonCorrelation(X, Y, log);
        return pearsonCorrelation;

    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A Pearson Correlation between two Parameters";
    }

    public String getExample() {
        return "<pearsonCorrelation id=\"r\" log=\"true\">\n"+
               "	<parameter idref=\"param1\"/>\n"+
               "    <parameter idref=\"param2\"/>\n"+
               "</pearsonCorrelation>\n";
    }

    public Class getReturnType() { return PearsonCorrelation.class; }

    public XMLSyntaxRule[] getSyntaxRules() { return rules; }

    private final XMLSyntaxRule[] rules = {
        // There should be two and only two Parameters (X & Y)
        new ElementRule(Parameter.class, 2, 2),
        // the optional log attribute has to be a Boolean
        AttributeRule.newBooleanRule(LOG, true),
    };

}

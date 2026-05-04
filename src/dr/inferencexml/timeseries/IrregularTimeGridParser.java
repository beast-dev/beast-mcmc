package dr.inferencexml.timeseries;

import dr.inference.timeseries.core.IrregularTimeGrid;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

/**
 * Parser for an explicit irregular time grid.
 */
public class IrregularTimeGridParser extends AbstractXMLObjectParser {

    public static final String PARSER_NAME = "irregularTimeGrid";
    public static final String TIMES = "times";

    @Override
    public String getParserName() {
        return PARSER_NAME;
    }

    @Override
    public Object parseXMLObject(final XMLObject xo) throws XMLParseException {
        final double[] times = xo.getDoubleArrayAttribute(TIMES);
        if (times.length < 1) {
            throw new XMLParseException("times must contain at least one entry");
        }
        for (int i = 1; i < times.length; ++i) {
            if (!(times[i] > times[i - 1])) {
                throw new XMLParseException("times must be strictly increasing");
            }
        }
        return new IrregularTimeGrid(times);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return RULES;
    }

    private static final XMLSyntaxRule[] RULES = new XMLSyntaxRule[] {
            AttributeRule.newDoubleArrayRule(TIMES)
    };

    @Override
    public String getParserDescription() {
        return "Defines an irregular time grid from explicit observation times.";
    }

    @Override
    public Class getReturnType() {
        return IrregularTimeGrid.class;
    }
}

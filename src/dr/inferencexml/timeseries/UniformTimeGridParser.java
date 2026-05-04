package dr.inferencexml.timeseries;

import dr.inference.timeseries.core.UniformTimeGrid;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

/**
 * Parser for a simple uniform time grid.
 */
public class UniformTimeGridParser extends AbstractXMLObjectParser {

    public static final String PARSER_NAME = "uniformTimeGrid";
    public static final String TIME_COUNT = "timeCount";
    public static final String START_TIME = "startTime";
    public static final String TIME_STEP = "timeStep";

    @Override
    public String getParserName() {
        return PARSER_NAME;
    }

    @Override
    public Object parseXMLObject(final XMLObject xo) throws XMLParseException {
        final int timeCount = xo.getIntegerAttribute(TIME_COUNT);
        final double startTime = xo.getAttribute(START_TIME, 0.0);
        final double timeStep = xo.getDoubleAttribute(TIME_STEP);

        if (timeCount < 1) {
            throw new XMLParseException("timeCount must be at least 1");
        }
        if (timeStep <= 0.0) {
            throw new XMLParseException("timeStep must be > 0.0");
        }

        return new UniformTimeGrid(timeCount, startTime, timeStep);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return RULES;
    }

    private static final XMLSyntaxRule[] RULES = new XMLSyntaxRule[] {
            AttributeRule.newIntegerRule(TIME_COUNT),
            AttributeRule.newDoubleRule(START_TIME, true),
            AttributeRule.newDoubleRule(TIME_STEP)
    };

    @Override
    public String getParserDescription() {
        return "Defines a simple regular time grid for time-series models.";
    }

    @Override
    public Class getReturnType() {
        return UniformTimeGrid.class;
    }
}

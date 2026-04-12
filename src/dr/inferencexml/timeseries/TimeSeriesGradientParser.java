package dr.inferencexml.timeseries;

import dr.inference.model.Parameter;
import dr.inference.timeseries.beast.TimeSeriesGradient;
import dr.inference.timeseries.likelihood.TimeSeriesLikelihood;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

/**
 * Parser for the BEAST-facing gradient adapter.
 */
public class TimeSeriesGradientParser extends AbstractXMLObjectParser {

    public static final String PARSER_NAME = "timeSeriesGradient";
    public static final String LIKELIHOOD = "likelihood";
    public static final String PARAMETER = "parameter";

    @Override
    public String getParserName() {
        return PARSER_NAME;
    }

    @Override
    public Object parseXMLObject(final XMLObject xo) throws XMLParseException {
        final TimeSeriesLikelihood likelihood =
                (TimeSeriesLikelihood) xo.getElementFirstChild(LIKELIHOOD);
        final Parameter parameter =
                (Parameter) xo.getElementFirstChild(PARAMETER);
        return new TimeSeriesGradient(likelihood, parameter);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return RULES;
    }

    private static final XMLSyntaxRule[] RULES = new XMLSyntaxRule[] {
            new ElementRule(LIKELIHOOD, new XMLSyntaxRule[] { new ElementRule(TimeSeriesLikelihood.class) }),
            new ElementRule(PARAMETER, new XMLSyntaxRule[] { new ElementRule(Parameter.class) })
    };

    @Override
    public String getParserDescription() {
        return "Creates a GradientWrtParameterProvider adapter for a time-series likelihood.";
    }

    @Override
    public Class getReturnType() {
        return TimeSeriesGradient.class;
    }
}

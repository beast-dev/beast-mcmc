package dr.inferencexml.timeseries;

import dr.inference.timeseries.likelihood.ParallelTimeSeriesLikelihood;
import dr.inference.timeseries.likelihood.TimeSeriesLikelihood;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

import java.util.ArrayList;
import java.util.List;

/**
 * Parser for independent time-series likelihoods sharing parameters.
 */
public class ParallelTimeSeriesLikelihoodParser extends AbstractXMLObjectParser {

    public static final String PARSER_NAME = "parallelTimeSeriesLikelihood";
    public static final String THREADS = "threads";

    @Override
    public String getParserName() {
        return PARSER_NAME;
    }

    @Override
    public Object parseXMLObject(final XMLObject xo) throws XMLParseException {
        final int threads = xo.getAttribute(THREADS, -1);
        if (threads < -1) {
            throw new XMLParseException(THREADS + " must be -1, 0, or a positive integer");
        }

        final List<TimeSeriesLikelihood> likelihoods = new ArrayList<TimeSeriesLikelihood>();
        for (int i = 0; i < xo.getChildCount(); ++i) {
            final Object child = xo.getChild(i);
            if (child instanceof TimeSeriesLikelihood) {
                if (likelihoods.contains(child)) {
                    throw new XMLParseException("The time-series likelihood '" +
                            ((TimeSeriesLikelihood) child).getId() + "' appears more than once.");
                }
                likelihoods.add((TimeSeriesLikelihood) child);
            } else {
                throw new XMLParseException("Only timeSeriesLikelihood children are allowed in "
                        + PARSER_NAME + ": found " + child);
            }
        }

        if (likelihoods.isEmpty()) {
            throw new XMLParseException(PARSER_NAME + " requires at least one timeSeriesLikelihood child");
        }

        final String id = xo.hasId() ? xo.getId() : PARSER_NAME;
        return new ParallelTimeSeriesLikelihood(id, threads, likelihoods);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return RULES;
    }

    private static final XMLSyntaxRule[] RULES = new XMLSyntaxRule[] {
            AttributeRule.newIntegerRule(THREADS, true),
            new ElementRule(TimeSeriesLikelihood.class, 1, Integer.MAX_VALUE)
    };

    @Override
    public String getParserDescription() {
        return "Aggregates independent time-series likelihoods and evaluates them in parallel.";
    }

    @Override
    public Class getReturnType() {
        return ParallelTimeSeriesLikelihood.class;
    }
}

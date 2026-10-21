package dr.inferencexml.timeseries;

import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.inference.timeseries.runtime.TimeSeriesBenchmark;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

/**
 * XML parser for timing already parsed time-series likelihood/gradient objects.
 */
public final class TimeSeriesBenchmarkParser extends AbstractXMLObjectParser {

    public static final String PARSER_NAME = "timeSeriesBenchmark";
    private static final String MODE = "mode";
    private static final String WARMUP = "warmup";
    private static final String N_WARMUP = "nWarmup";
    private static final String TIMED_ITERATIONS = "timedIterations";
    private static final String N_TIMED = "nTimed";
    private static final String ITERATION_COUNT = "iterationCount";
    private static final String PERTURB_EACH_ITERATION = "perturbEachIteration";
    private static final String PERTURBATION = "perturbation";
    private static final String PERTURB_INDEX = "perturbIndex";
    private static final String DESCRIPTION = "description";

    @Override
    public String getParserName() {
        return PARSER_NAME;
    }

    @Override
    public Object parseXMLObject(final XMLObject xo) throws XMLParseException {
        final TimeSeriesBenchmark.Mode mode;
        try {
            mode = TimeSeriesBenchmark.Mode.parse(xo.getAttribute(MODE, "gradient"));
        } catch (IllegalArgumentException e) {
            throw new XMLParseException(e.getMessage());
        }

        final Likelihood likelihood = (Likelihood) xo.getChild(Likelihood.class);
        final GradientWrtParameterProvider gradient =
                (GradientWrtParameterProvider) xo.getChild(GradientWrtParameterProvider.class);
        final boolean perturbEachIteration = xo.getAttribute(PERTURB_EACH_ITERATION, true);
        final Parameter perturbParameter = (Parameter) xo.getChild(Parameter.class);

        if (likelihood == null) {
            throw new XMLParseException(PARSER_NAME + " requires a likelihood child");
        }
        if (mode == TimeSeriesBenchmark.Mode.GRADIENT && gradient == null) {
            throw new XMLParseException(PARSER_NAME + " requires a gradient child in gradient mode");
        }
        if (perturbEachIteration && perturbParameter == null) {
            throw new XMLParseException(PARSER_NAME + " requires a parameter child when perturbEachIteration=true");
        }

        return new TimeSeriesBenchmark(
                xo.hasId() ? xo.getId() : PARSER_NAME,
                mode,
                likelihood,
                gradient,
                getAliasedIntegerAttribute(xo, WARMUP, N_WARMUP, 20),
                getTimedIterations(xo),
                perturbEachIteration,
                perturbParameter,
                xo.getAttribute(PERTURB_INDEX, 0),
                xo.getAttribute(PERTURBATION, 1.0e-9),
                xo.getAttribute(DESCRIPTION, ""));
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return RULES;
    }

    private static final XMLSyntaxRule[] RULES = new XMLSyntaxRule[] {
            AttributeRule.newStringRule(MODE, true),
            AttributeRule.newIntegerRule(WARMUP, true),
            AttributeRule.newIntegerRule(N_WARMUP, true),
            AttributeRule.newIntegerRule(TIMED_ITERATIONS, true),
            AttributeRule.newIntegerRule(N_TIMED, true),
            AttributeRule.newIntegerRule(ITERATION_COUNT, true),
            AttributeRule.newBooleanRule(PERTURB_EACH_ITERATION, true),
            AttributeRule.newDoubleRule(PERTURBATION, true),
            AttributeRule.newIntegerRule(PERTURB_INDEX, true),
            AttributeRule.newStringRule(DESCRIPTION, true),
            new ElementRule(Likelihood.class),
            new ElementRule(GradientWrtParameterProvider.class, true),
            new ElementRule(Parameter.class, true)
    };

    @Override
    public String getParserDescription() {
        return "Times already parsed time-series likelihood and gradient objects.";
    }

    @Override
    public Class getReturnType() {
        return TimeSeriesBenchmark.class;
    }

    private static int getTimedIterations(final XMLObject xo) throws XMLParseException {
        if (xo.hasAttribute(TIMED_ITERATIONS)) {
            return xo.getIntegerAttribute(TIMED_ITERATIONS);
        }
        if (xo.hasAttribute(N_TIMED)) {
            return xo.getIntegerAttribute(N_TIMED);
        }
        if (xo.hasAttribute(ITERATION_COUNT)) {
            return xo.getIntegerAttribute(ITERATION_COUNT);
        }
        return 100;
    }

    private static int getAliasedIntegerAttribute(final XMLObject xo,
                                                  final String primary,
                                                  final String alias,
                                                  final int defaultValue) throws XMLParseException {
        if (xo.hasAttribute(primary)) {
            return xo.getIntegerAttribute(primary);
        }
        if (xo.hasAttribute(alias)) {
            return xo.getIntegerAttribute(alias);
        }
        return defaultValue;
    }
}

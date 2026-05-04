package dr.inferencexml.timeseries;

import dr.inference.timeseries.core.TimeSeriesModel;
import dr.inference.timeseries.gaussian.GaussianObservationModel;
import dr.inference.timeseries.gaussian.OUTimeSeriesProcessAdapter;
import dr.inference.timeseries.engine.gaussian.GaussianForwardComputationMode;
import dr.inference.timeseries.likelihood.GaussianGradientComputationMode;
import dr.inference.timeseries.likelihood.GaussianSmootherComputationMode;
import dr.inference.timeseries.likelihood.GaussianTimeSeriesLikelihoodFactory;
import dr.inference.timeseries.likelihood.TimeSeriesLikelihood;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

/**
 * Parser constructing a time-series likelihood for the OU plus Gaussian observation path.
 */
public class TimeSeriesLikelihoodParser extends AbstractXMLObjectParser {

    public static final String PARSER_NAME = "timeSeriesLikelihood";
    public static final String MODEL = "model";
    public static final String FORWARD_MODE = "forwardMode";
    public static final String SMOOTHER_MODE = "smootherMode";
    public static final String GRADIENT_MODE = "gradientMode";
    public static final String DEBUG_MODES = "debugModes";

    @Override
    public String getParserName() {
        return PARSER_NAME;
    }

    @Override
    public Object parseXMLObject(final XMLObject xo) throws XMLParseException {
        final TimeSeriesModel model = (TimeSeriesModel) xo.getElementFirstChild(MODEL);

        if (!(model.getLatentProcessModel() instanceof OUTimeSeriesProcessAdapter)) {
            throw new XMLParseException("This initial parser currently supports OUTimeSeriesProcessAdapter only.");
        }
        if (!(model.getObservationModel() instanceof GaussianObservationModel)) {
            throw new XMLParseException("This initial parser currently supports GaussianObservationModel only.");
        }

        final String id = xo.hasId() ? xo.getId() : PARSER_NAME;
        final GaussianForwardComputationMode forwardMode = parseForwardMode(xo);
        final GaussianSmootherComputationMode smootherMode = parseSmootherMode(xo);
        final GaussianGradientComputationMode gradientMode = parseGradientMode(xo);
        validateProductionModes(xo, forwardMode, smootherMode, gradientMode);

        return GaussianTimeSeriesLikelihoodFactory.create(
                id,
                model,
                forwardMode,
                smootherMode,
                gradientMode);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return RULES;
    }

    private static final XMLSyntaxRule[] RULES = new XMLSyntaxRule[] {
            AttributeRule.newStringRule(FORWARD_MODE, true,
                    "Timeseries forward mode. Defaults to canonical; expectation is debug-only."),
            AttributeRule.newStringRule(SMOOTHER_MODE, true,
                    "Timeseries smoother mode. Defaults to canonical; expectation is debug-only."),
            AttributeRule.newStringRule(GRADIENT_MODE, true,
                    "Timeseries gradient mode. Defaults to canonicalAnalytical; other values are debug-only."),
            AttributeRule.newBooleanRule(DEBUG_MODES, true,
                    "Must be true to use expectation or disabled debug modes."),
            new ElementRule(MODEL, new XMLSyntaxRule[] { new ElementRule(TimeSeriesModel.class) })
    };

    @Override
    public String getParserDescription() {
        return "Constructs a time-series likelihood for the initial OU plus Gaussian observation path.";
    }

    @Override
    public Class getReturnType() {
        return TimeSeriesLikelihood.class;
    }

    private static GaussianForwardComputationMode parseForwardMode(final XMLObject xo)
            throws XMLParseException {
        final String value = xo.getAttribute(FORWARD_MODE, "canonical");
        if ("canonical".equalsIgnoreCase(value)) {
            return GaussianForwardComputationMode.CANONICAL;
        }
        if ("expectation".equalsIgnoreCase(value)) {
            return GaussianForwardComputationMode.EXPECTATION;
        }
        throw new XMLParseException("Unknown " + FORWARD_MODE + ": " + value);
    }

    private static GaussianSmootherComputationMode parseSmootherMode(final XMLObject xo)
            throws XMLParseException {
        final String value = xo.getAttribute(SMOOTHER_MODE, "canonical");
        if ("canonical".equalsIgnoreCase(value)) {
            return GaussianSmootherComputationMode.CANONICAL;
        }
        if ("expectation".equalsIgnoreCase(value)) {
            return GaussianSmootherComputationMode.EXPECTATION;
        }
        throw new XMLParseException("Unknown " + SMOOTHER_MODE + ": " + value);
    }

    private static GaussianGradientComputationMode parseGradientMode(final XMLObject xo)
            throws XMLParseException {
        final String value = xo.getAttribute(GRADIENT_MODE, "canonicalAnalytical");
        if ("canonicalAnalytical".equalsIgnoreCase(value)
                || "canonical".equalsIgnoreCase(value)) {
            return GaussianGradientComputationMode.CANONICAL_ANALYTICAL;
        }
        if ("expectationAnalytical".equalsIgnoreCase(value)
                || "expectation".equalsIgnoreCase(value)) {
            return GaussianGradientComputationMode.EXPECTATION_ANALYTICAL;
        }
        if ("disabled".equalsIgnoreCase(value)
                || "none".equalsIgnoreCase(value)) {
            return GaussianGradientComputationMode.DISABLED;
        }
        throw new XMLParseException("Unknown " + GRADIENT_MODE + ": " + value);
    }

    private static void validateProductionModes(final XMLObject xo,
                                                final GaussianForwardComputationMode forwardMode,
                                                final GaussianSmootherComputationMode smootherMode,
                                                final GaussianGradientComputationMode gradientMode)
            throws XMLParseException {
        final boolean canonicalProduction =
                forwardMode == GaussianForwardComputationMode.CANONICAL
                        && smootherMode == GaussianSmootherComputationMode.CANONICAL
                        && gradientMode == GaussianGradientComputationMode.CANONICAL_ANALYTICAL;
        if (canonicalProduction || xo.getAttribute(DEBUG_MODES, false)) {
            return;
        }
        throw new XMLParseException("Timeseries likelihood uses canonical mode by default. "
                + "Set " + DEBUG_MODES + "=\"true\" to use expectation or disabled debug modes.");
    }
}

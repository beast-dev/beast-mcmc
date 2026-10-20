package dr.inferencexml.timeseries;

import dr.inference.model.MatrixParameter;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Parameter;
import dr.evomodel.continuous.ou.OUProcessModel;
import dr.inference.timeseries.model.gaussian.OUTimeSeriesProcessAdapter;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.ElementRule;
import dr.xml.StringAttributeRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

/**
 * Parser for the OU process model scaffold.
 */
public class OUProcessModelParser extends AbstractXMLObjectParser {

    public static final String PARSER_NAME = "ouProcessModel";
    public static final String STATE_DIMENSION = "stateDimension";
    public static final String DRIFT_MATRIX = "driftMatrix";
    public static final String DIFFUSION_MATRIX = "diffusionMatrix";
    public static final String STATIONARY_MEAN = "stationaryMean";
    public static final String INITIAL_COVARIANCE = "initialCovariance";
    public static final String DIFFUSION_FLOOR = "diffusionFloor";

    @Override
    public String getParserName() {
        return PARSER_NAME;
    }

    @Override
    public Object parseXMLObject(final XMLObject xo) throws XMLParseException {
        final int stateDimension = xo.getIntegerAttribute(STATE_DIMENSION);
        final MatrixParameterInterface driftMatrix =
                (MatrixParameterInterface) xo.getElementFirstChild(DRIFT_MATRIX);
        final MatrixParameter diffusionMatrix = (MatrixParameter) xo.getElementFirstChild(DIFFUSION_MATRIX);
        final Parameter stationaryMean = (Parameter) xo.getElementFirstChild(STATIONARY_MEAN);
        final MatrixParameter initialCovariance = (MatrixParameter) xo.getElementFirstChild(INITIAL_COVARIANCE);

        OUSelectionChartParserHelper.validateSelectionChart(xo, driftMatrix, PARSER_NAME);

        final double diffusionFloor = xo.getAttribute(DIFFUSION_FLOOR, 0.0);
        if (!Double.isFinite(diffusionFloor) || diffusionFloor < 0.0) {
            throw new XMLParseException(DIFFUSION_FLOOR + " must be a finite non-negative number, found "
                    + diffusionFloor);
        }

        final String id = xo.hasId() ? xo.getId() : PARSER_NAME;
        return new OUTimeSeriesProcessAdapter(
                new OUProcessModel(id, stateDimension, driftMatrix, diffusionMatrix, stationaryMean,
                        initialCovariance, diffusionFloor));
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return RULES;
    }

    private static final XMLSyntaxRule[] RULES = new XMLSyntaxRule[] {
            AttributeRule.newIntegerRule(STATE_DIMENSION),
            AttributeRule.newDoubleRule(DIFFUSION_FLOOR, true,
                    "Constant added to every diagonal entry of the diffusion matrix (Q + floor I) in "
                            + "per-time-unit variance; default 0 leaves the model unchanged."),
            new StringAttributeRule(OUSelectionChartParserHelper.SELECTION_CHART,
                    "Selection-matrix chart for OU models. Orthogonal block is the default; dense must be explicit.",
                    OUSelectionChartParserHelper.ALLOWED_SELECTION_CHARTS, true),
            new ElementRule(DRIFT_MATRIX, new XMLSyntaxRule[] { new ElementRule(MatrixParameterInterface.class) }),
            new ElementRule(DIFFUSION_MATRIX, new XMLSyntaxRule[] { new ElementRule(MatrixParameter.class) }),
            new ElementRule(STATIONARY_MEAN, new XMLSyntaxRule[] { new ElementRule(Parameter.class) }),
            new ElementRule(INITIAL_COVARIANCE, new XMLSyntaxRule[] { new ElementRule(MatrixParameter.class) })
    };

    @Override
    public String getParserDescription() {
        return "Defines a multivariate OU latent process adapter for Gaussian time-series inference.";
    }

    @Override
    public Class getReturnType() {
        return OUTimeSeriesProcessAdapter.class;
    }
}

package dr.inferencexml.timeseries;

import dr.inference.model.MatrixParameter;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Parameter;
import dr.evomodel.continuous.ou.OUProcessModel;
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

        final String id = xo.hasId() ? xo.getId() : PARSER_NAME;
        return new OUProcessModel(id, stateDimension, driftMatrix, diffusionMatrix, stationaryMean, initialCovariance);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return RULES;
    }

    private static final XMLSyntaxRule[] RULES = new XMLSyntaxRule[] {
            AttributeRule.newIntegerRule(STATE_DIMENSION),
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
        return "Defines a multivariate OU latent process model with a Gaussian transition representation.";
    }

    @Override
    public Class getReturnType() {
        return OUProcessModel.class;
    }
}

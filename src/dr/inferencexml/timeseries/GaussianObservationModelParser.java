package dr.inferencexml.timeseries;

import dr.inference.model.MatrixParameter;
import dr.inference.timeseries.gaussian.GaussianObservationModel;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

/**
 * Parser for the Gaussian observation model scaffold.
 */
public class GaussianObservationModelParser extends AbstractXMLObjectParser {

    public static final String PARSER_NAME = "gaussianObservationModel";
    public static final String OBSERVATION_DIMENSION = "observationDimension";
    public static final String DESIGN_MATRIX = "designMatrix";
    public static final String NOISE_COVARIANCE = "noiseCovariance";
    public static final String OBSERVATIONS = "observations";

    @Override
    public String getParserName() {
        return PARSER_NAME;
    }

    @Override
    public Object parseXMLObject(final XMLObject xo) throws XMLParseException {
        final int observationDimension = xo.getIntegerAttribute(OBSERVATION_DIMENSION);
        final MatrixParameter designMatrix = (MatrixParameter) xo.getElementFirstChild(DESIGN_MATRIX);
        final MatrixParameter noiseCovariance = (MatrixParameter) xo.getElementFirstChild(NOISE_COVARIANCE);
        final MatrixParameter observations = (MatrixParameter) xo.getElementFirstChild(OBSERVATIONS);

        final String id = xo.hasId() ? xo.getId() : PARSER_NAME;
        return new GaussianObservationModel(id, observationDimension, designMatrix, noiseCovariance, observations);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return RULES;
    }

    private static final XMLSyntaxRule[] RULES = new XMLSyntaxRule[] {
            AttributeRule.newIntegerRule(OBSERVATION_DIMENSION),
            new ElementRule(DESIGN_MATRIX, new XMLSyntaxRule[] { new ElementRule(MatrixParameter.class) }),
            new ElementRule(NOISE_COVARIANCE, new XMLSyntaxRule[] { new ElementRule(MatrixParameter.class) }),
            new ElementRule(OBSERVATIONS, new XMLSyntaxRule[] { new ElementRule(MatrixParameter.class) })
    };

    @Override
    public String getParserDescription() {
        return "Defines a Gaussian observation model for a latent Gaussian time-series process.";
    }

    @Override
    public Class getReturnType() {
        return GaussianObservationModel.class;
    }
}

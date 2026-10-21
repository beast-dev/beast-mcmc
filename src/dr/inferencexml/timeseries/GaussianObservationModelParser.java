package dr.inferencexml.timeseries;

import dr.inference.model.MatrixParameter;
import dr.inference.timeseries.model.gaussian.LinearGaussianObservationModel;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

/**
 * Parser for the linear Gaussian observation model.
 *
 * <p>Matrix shapes are validated here, at parse time, so that a wrongly shaped matrix (most
 * often an observation matrix written with one {@code <parameter>} per observed coordinate
 * instead of one per time point) fails with a message naming the matrix rather than with an
 * index error during the first likelihood evaluation.  The design matrix column count is
 * checked against the latent state dimension by {@link BasicTimeSeriesModelParser}, where
 * the latent process is known.
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
        if (observationDimension < 1) {
            throw new XMLParseException(OBSERVATION_DIMENSION + " must be at least 1; found " + observationDimension);
        }
        final MatrixParameter designMatrix = (MatrixParameter) xo.getElementFirstChild(DESIGN_MATRIX);
        final MatrixParameter noiseCovariance = (MatrixParameter) xo.getElementFirstChild(NOISE_COVARIANCE);
        final MatrixParameter observations = (MatrixParameter) xo.getElementFirstChild(OBSERVATIONS);

        final String id = xo.hasId() ? xo.getId() : PARSER_NAME;
        try {
            return new LinearGaussianObservationModel(id, observationDimension, designMatrix, noiseCovariance, observations);
        } catch (final IllegalArgumentException e) {
            // The constructor checks the matrix shapes; surface its message as a parser error.
            throw new XMLParseException(e.getMessage());
        }
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
        return "Defines a linear Gaussian observation model y_t = H x_t + e_t for a latent Gaussian"
                + " time-series process. designMatrix H is observationDimension x stateDimension,"
                + " noiseCovariance R is observationDimension x observationDimension, and observations"
                + " is observationDimension x timeCount with one matrixParameter column (one <parameter>"
                + " child) per time point; NaN cells mark missing values.";
    }

    @Override
    public Class getReturnType() {
        return LinearGaussianObservationModel.class;
    }
}

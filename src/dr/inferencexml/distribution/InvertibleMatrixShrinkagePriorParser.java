package dr.inferencexml.distribution;

import dr.inference.distribution.InvertibleMatrixShrinkagePrior;
import dr.inference.model.MatrixParameter;
import dr.inference.model.Parameter;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

/**
 * Parses {@link InvertibleMatrixShrinkagePrior}.
 */
public class InvertibleMatrixShrinkagePriorParser extends AbstractXMLObjectParser {

    public static final String INVERTIBLE_MATRIX_SHRINKAGE_PRIOR = "invertibleMatrixShrinkagePrior";
    private static final String MATRIX = "matrix";
    private static final String LAMBDA = "lambda";
    private static final String SCALE_FROBENIUS_BY_DIMENSION = "scaleFrobeniusByDimension";

    @Override
    public String getParserName() {
        return INVERTIBLE_MATRIX_SHRINKAGE_PRIOR;
    }

    @Override
    public Object parseXMLObject(final XMLObject xo) throws XMLParseException {

        final MatrixParameter matrix = (MatrixParameter)
                xo.getElementFirstChild(MATRIX);
        final Parameter lambda = (Parameter) xo.getElementFirstChild(LAMBDA);
        final boolean scaleFrobeniusByDimension =
                xo.getAttribute(SCALE_FROBENIUS_BY_DIMENSION, false);
        final String name = xo.hasId() ? xo.getId() : getParserName();

        return new InvertibleMatrixShrinkagePrior(name, matrix, lambda, scaleFrobeniusByDimension);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newBooleanRule(SCALE_FROBENIUS_BY_DIMENSION, true,
                    "If true, use log p(R) = C + lambda * log|det R| "
                            + "- (lambda * dim(R) / 2) ||R||_F^2. Default false uses "
                            + "- (lambda / 2) ||R||_F^2."),
            new ElementRule(MATRIX, new XMLSyntaxRule[]{
                    new ElementRule(MatrixParameter.class)
            }),
            new ElementRule(LAMBDA, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }),
    };

    @Override
    public String getParserDescription() {
        return "A prior on a square matrix R that repels R from singularity and bounds its "
                + "overall scale: log p(R) = C + lambda * log|det R| "
                + "- (lambda * c / 2) ||R||_F^2, where c = 1 by default and "
                + "c = dim(R) when scaleFrobeniusByDimension is true.";
    }

    @Override
    public Class getReturnType() {
        return InvertibleMatrixShrinkagePrior.class;
    }
}

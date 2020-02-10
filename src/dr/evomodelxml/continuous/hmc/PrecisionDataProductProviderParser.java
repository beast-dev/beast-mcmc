package dr.evomodelxml.continuous.hmc;

import dr.inference.hmc.PrecisionMatrixVectorProductProvider;
import dr.inference.model.MatrixParameter;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Parameter;
import dr.inferencexml.model.MaskedParameterParser;
import dr.xml.*;

/**
 * @author Zhenyu Zhang
 */

public class PrecisionDataProductProviderParser extends AbstractXMLObjectParser {
    private static final String PRODUCT_PROVIDER = "precisionVectorProduct";
    private static final String MASKING = MaskedParameterParser.MASKING;
    private static final String TIME_GUESS = "roughTravelTimeGuess";
    private static final String THREAD_COUNT = "threadCount";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        double roughTimeGuess = xo.getAttribute(TIME_GUESS, 1.0);
        MatrixParameter matrix = (MatrixParameter) xo.getChild(MatrixParameterInterface.class);
        return new PrecisionMatrixVectorProductProvider.Generic(matrix, roughTimeGuess);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newIntegerRule(THREAD_COUNT, true),
            AttributeRule.newDoubleRule(TIME_GUESS, true),
            new ElementRule(MASKING,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class)
                    }, true),
    };

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return PrecisionMatrixVectorProductProvider.class;
    }

    @Override
    public String getParserName() {
        return PRODUCT_PROVIDER;
    }
}

package dr.inferencexml.distribution;

import dr.inference.distribution.RowDimensionMultinomialPrior;
import dr.inference.model.AdaptableSizeFastMatrixParameter;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * @author Max Tolkoff
 */
public class RowDimensionMultinomialPriorParser extends AbstractXMLObjectParser {
    public static final String ROW_DIMENSION_MULTINOMIAL_PRIOR = "rowDimensionMultinomialPrior";
    public static final String PROBABILITIES = "probabilities";
    public static final String TRANSPOSE = "transpose";


    @Override
    public String getParserName() {
        return ROW_DIMENSION_MULTINOMIAL_PRIOR;
    }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        AdaptableSizeFastMatrixParameter data = (AdaptableSizeFastMatrixParameter) xo.getChild(AdaptableSizeFastMatrixParameter.class);
        Parameter probabilities = (Parameter) xo.getChild(PROBABILITIES).getChild(Parameter.class);
        boolean transpose = false;
        if(xo.hasAttribute(TRANSPOSE))
            transpose = xo.getBooleanAttribute(TRANSPOSE);
        String id = xo.getId();


        return new RowDimensionMultinomialPrior(id, data, probabilities, transpose);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(PROBABILITIES,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class)
                    }
            ),
            new ElementRule(AdaptableSizeFastMatrixParameter.class),
            AttributeRule.newBooleanRule(TRANSPOSE, true),
    };


    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return RowDimensionMultinomialPrior.class;
    }
}

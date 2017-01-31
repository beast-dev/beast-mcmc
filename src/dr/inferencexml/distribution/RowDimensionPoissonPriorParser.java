package dr.inferencexml.distribution;


import dr.inference.distribution.RowDimensionPoissonPrior;
import dr.inference.model.AdaptableSizeFastMatrixParameter;
import dr.xml.*;

/**
 * Created by maxryandolinskytolkoff on 7/20/16.
 */
public class RowDimensionPoissonPriorParser extends AbstractXMLObjectParser{
    public static final String ROW_DIMENSION_POISSON_PRIOR = "rowDimensionPoissonPrior";
    public static final String UNTRUNCATED_MEAN = "untruncatedMean";
    public static final String TRANSPOSE = "transpose";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        AdaptableSizeFastMatrixParameter parameter = (AdaptableSizeFastMatrixParameter) xo.getChild(AdaptableSizeFastMatrixParameter.class);
        double untruncatedMean = xo.getDoubleAttribute(UNTRUNCATED_MEAN);
        boolean transpose = false;
        if(xo.hasAttribute(TRANSPOSE))
            transpose = xo.getBooleanAttribute(TRANSPOSE);

        return new RowDimensionPoissonPrior(untruncatedMean, parameter, transpose);
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            AttributeRule.newDoubleRule(UNTRUNCATED_MEAN, true),
            AttributeRule.newBooleanRule(TRANSPOSE, true),
            new ElementRule(AdaptableSizeFastMatrixParameter.class)
    };

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return RowDimensionPoissonPrior.class;
    }

    @Override
    public String getParserName() {
        return ROW_DIMENSION_POISSON_PRIOR;
    }
}

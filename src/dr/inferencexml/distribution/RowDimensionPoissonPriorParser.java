package dr.inferencexml.distribution;


import dr.inference.distribution.DeterminentalPointProcessPrior;
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
        AdaptableSizeFastMatrixParameter parameter;
        DeterminentalPointProcessPrior DPP;
        if(xo.getChild(AdaptableSizeFastMatrixParameter.class) != null)
           parameter= (AdaptableSizeFastMatrixParameter) xo.getChild(AdaptableSizeFastMatrixParameter.class);
        else
            parameter= null;
        if(xo.getChild(DeterminentalPointProcessPrior.class) != null)
            DPP = (DeterminentalPointProcessPrior) xo.getChild(DeterminentalPointProcessPrior.class);
        else
            DPP = null;
        double untruncatedMean = xo.getDoubleAttribute(UNTRUNCATED_MEAN);
        boolean transpose = false;
        if(xo.hasAttribute(TRANSPOSE))
            transpose = xo.getBooleanAttribute(TRANSPOSE);
        String id = xo.getId();

        return new RowDimensionPoissonPrior(id, untruncatedMean, parameter, DPP, transpose);
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            AttributeRule.newDoubleRule(UNTRUNCATED_MEAN, true),
            AttributeRule.newBooleanRule(TRANSPOSE, true),
            new OrRule(
            new ElementRule(AdaptableSizeFastMatrixParameter.class),
                    new ElementRule(DeterminentalPointProcessPrior.class)
            )
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

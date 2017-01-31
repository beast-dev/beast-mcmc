package dr.inferencexml.operators;

import dr.evomodel.continuous.OrderedLatentLiabilityLikelihood;
import dr.inference.model.LatentFactorModel;
import dr.inference.operators.LatentFactorLiabilityGibbsOperator;
import dr.xml.*;

/**
 * Created by Max on 9/1/16.
 */
public class LatentFactorLiabilityGibbsOperatorParser extends AbstractXMLObjectParser{
    public final static String LATENT_FACTOR_LIABILITY_GIBBS_OPERATOR = "latentFactorLiabilityGibbsOperator";
    public final static String WEIGHT = "weight";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        double weight = xo.getDoubleAttribute(WEIGHT);
        LatentFactorModel lfm = (LatentFactorModel) xo.getChild(LatentFactorModel.class);
        OrderedLatentLiabilityLikelihood liabilityLikelihood = (OrderedLatentLiabilityLikelihood) xo.getChild(OrderedLatentLiabilityLikelihood.class);

        return new LatentFactorLiabilityGibbsOperator(weight, lfm, liabilityLikelihood);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            AttributeRule.newDoubleRule(WEIGHT),
            new ElementRule(LatentFactorModel.class),
            new ElementRule(OrderedLatentLiabilityLikelihood.class)
    };

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return LatentFactorLiabilityGibbsOperator.class;
    }

    @Override
    public String getParserName() {
        return LATENT_FACTOR_LIABILITY_GIBBS_OPERATOR;
    }
}

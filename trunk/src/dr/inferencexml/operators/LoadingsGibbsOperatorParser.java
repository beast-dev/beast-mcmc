package dr.inferencexml.operators;

import dr.evomodel.continuous.LatentFactorModel;
import dr.inference.distribution.MultivariateDistributionLikelihood;
import dr.inference.operators.LoadingsGibbsOperator;
import dr.xml.*;
import dr.inference.distribution.DistributionLikelihood;

/**
 * Created with IntelliJ IDEA.
 * User: max
 * Date: 5/23/14
 * Time: 1:49 PM
 * To change this template use File | Settings | File Templates.
 */
public class LoadingsGibbsOperatorParser extends AbstractXMLObjectParser {
    public static final String LOADINGS_GIBBS_OPERATOR="loadingsGibbsOperator";
    public static final String WEIGHT="weight";


    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        String weightTemp= (String) xo.getAttribute(WEIGHT);
        Double weight=Double.parseDouble(weightTemp);
        LatentFactorModel LFM =(LatentFactorModel) xo.getChild(LatentFactorModel.class);
        DistributionLikelihood prior= (DistributionLikelihood) xo.getChild(DistributionLikelihood.class);


        return new LoadingsGibbsOperator(LFM, prior, weight);  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(LatentFactorModel.class),
            new ElementRule(DistributionLikelihood.class),
//            new ElementRule(CompoundParameter.class),
            AttributeRule.newDoubleRule(WEIGHT),
    };

    @Override
    public String getParserDescription() {
        return "Gibbs sampler for the loadings matrix of a latent factor model";  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public Class getReturnType() {
        return LoadingsGibbsOperator.class;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public String getParserName() {
        return LOADINGS_GIBBS_OPERATOR;  //To change body of implemented methods use File | Settings | File Templates.
    }
}

package dr.inferencexml.operators;

import dr.evomodel.continuous.LatentFactorModel;
import dr.inference.distribution.DistributionLikelihood;
import dr.inference.operators.LatentFactorModelPrecisionGibbsOperator;
import dr.inference.operators.LoadingsGibbsOperator;
import dr.xml.*;

/**
 * Created by max on 6/12/14.
 */
public class LatentFactorModelPrecisionGibbsOperatorParser extends AbstractXMLObjectParser {
    public final String LATENT_FACTOR_MODEL_PRECISION_OPERATOR="latentFactorModelPrecisionOperator";
    public final String WEIGHT="weight";
    public final String RANDOM_SCAN="randomScan";


    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        String weightTemp= (String) xo.getAttribute(WEIGHT);
        double weight=Double.parseDouble(weightTemp);
        LatentFactorModel LFM =(LatentFactorModel) xo.getChild(LatentFactorModel.class);
        DistributionLikelihood prior= (DistributionLikelihood) xo.getChild(DistributionLikelihood.class);
        boolean randomScan=xo.getAttribute(RANDOM_SCAN, true);

        return new LatentFactorModelPrecisionGibbsOperator(LFM, prior, weight, randomScan);


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
        return "Gibbs sampler for the precision of a factor analysis model";
    }

    @Override
    public Class getReturnType() {
        return null;
    }

    @Override
    public String getParserName() {
        return LATENT_FACTOR_MODEL_PRECISION_OPERATOR;
    }
}

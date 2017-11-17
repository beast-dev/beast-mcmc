package dr.inferencexml.operators;

import dr.inference.distribution.DistributionLikelihood;
import dr.inference.distribution.IndependentInverseGammaDistributionModel;
import dr.inference.distribution.ShrinkageGibbsOperator;
import dr.inference.model.Parameter;
import dr.inference.operators.ShrinkageAugmentedGibbsOperator;
import dr.xml.*;

public class ShrinkageAugmentedGibbsOperatorParser extends AbstractXMLObjectParser{
    public final static String LOCAL = "local";
    public final static String GLOBAL = "global";
    public final static String SHRINKAGE_AUGMENTED_GIBBS_OPERATOR = "shrinkageAugmentedGibbsOperator";
    public final static String WEIGHT = "weight";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        double weight = xo.getDoubleAttribute(WEIGHT);
        XMLObject local = xo.getChild(LOCAL);
        XMLObject global = xo.getChild(GLOBAL);

        IndependentInverseGammaDistributionModel localPrior = (IndependentInverseGammaDistributionModel) local.getChild(IndependentInverseGammaDistributionModel.class);
        DistributionLikelihood localAugmentedPrior = (DistributionLikelihood) local.getChild(DistributionLikelihood.class);

        IndependentInverseGammaDistributionModel globalPrior = (IndependentInverseGammaDistributionModel) global.getChild(IndependentInverseGammaDistributionModel.class);
        DistributionLikelihood globalAugmentedPrior = (DistributionLikelihood) global.getChild(DistributionLikelihood.class);



        return new ShrinkageAugmentedGibbsOperator(weight, localAugmentedPrior, globalAugmentedPrior, localPrior, globalPrior);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            AttributeRule.newDoubleRule(WEIGHT),
            new ElementRule(LOCAL, new XMLSyntaxRule[]{
                    new ElementRule(IndependentInverseGammaDistributionModel.class),
                    new ElementRule(DistributionLikelihood.class)
            }),
            new ElementRule(GLOBAL, new XMLSyntaxRule[]{
                    new ElementRule(IndependentInverseGammaDistributionModel.class),
                    new ElementRule(DistributionLikelihood.class)
            }),
    };

    @Override
    public String getParserDescription() {
        return "Gibbs sampler for augmented portions of a shrinkage prior";
    }

    @Override
    public Class getReturnType() {
        return ShrinkageGibbsOperator.class;
    }

    @Override
    public String getParserName() {
        return SHRINKAGE_AUGMENTED_GIBBS_OPERATOR;
    }
}

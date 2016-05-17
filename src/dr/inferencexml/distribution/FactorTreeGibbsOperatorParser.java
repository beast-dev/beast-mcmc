package dr.inferencexml.distribution;

import dr.evomodel.continuous.FullyConjugateMultivariateTraitLikelihood;
import dr.inference.distribution.DistributionLikelihood;
import dr.inference.distribution.MomentDistributionModel;
import dr.inference.model.LatentFactorModel;
import dr.inference.operators.FactorTreeGibbsOperator;
import dr.xml.*;

/**
 * Created by max on 5/16/16.
 */
public class FactorTreeGibbsOperatorParser extends AbstractXMLObjectParser {
    public static final String FACTOR_TREE_GIBBS_OPERATOR_PARSER = "factorTreeGibbsOperator";
    public static final String WEIGHT = "weight";


    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        double weight = (Double) xo.getAttribute(WEIGHT);
        LatentFactorModel lfm= (LatentFactorModel) xo.getChild(LatentFactorModel.class);
        FullyConjugateMultivariateTraitLikelihood tree = (FullyConjugateMultivariateTraitLikelihood) xo.getChild(FullyConjugateMultivariateTraitLikelihood.class);

        return new FactorTreeGibbsOperator(weight, lfm, tree);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[0];
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(LatentFactorModel.class),
            new ElementRule(FullyConjugateMultivariateTraitLikelihood.class),
            AttributeRule.newDoubleRule(WEIGHT),
    };


    @Override
    public String getParserDescription() {
        return "Gibbs sample a factor column given a FullyConjugateMultivariateTree";
    }

    @Override
    public Class getReturnType() {
        return FactorTreeGibbsOperator.class;
    }

    @Override
    public String getParserName() {
        return FACTOR_TREE_GIBBS_OPERATOR_PARSER;
    }
}

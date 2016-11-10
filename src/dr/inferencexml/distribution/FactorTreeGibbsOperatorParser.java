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
    public static final String RANDOM_SCAN = "randomScan";
    public static final String WORKING_PRIOR = "workingPrior";


    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        String weightTemp = (String) xo.getAttribute(WEIGHT);
        double weight = Double.parseDouble(weightTemp);
        LatentFactorModel lfm= (LatentFactorModel) xo.getChild(LatentFactorModel.class);
        FullyConjugateMultivariateTraitLikelihood tree = (FullyConjugateMultivariateTraitLikelihood) xo.getChild(FullyConjugateMultivariateTraitLikelihood.class);
        FullyConjugateMultivariateTraitLikelihood workingTree = null;
        if(xo.getChild(WORKING_PRIOR) != null){
            System.out.println("happy");
            workingTree = (FullyConjugateMultivariateTraitLikelihood) xo.getChild(WORKING_PRIOR).getChild(FullyConjugateMultivariateTraitLikelihood.class);
        }
        boolean randomScan = xo.getAttribute(RANDOM_SCAN, true);

        return new FactorTreeGibbsOperator(weight, lfm, tree, randomScan);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(LatentFactorModel.class),
            new ElementRule(FullyConjugateMultivariateTraitLikelihood.class),
            AttributeRule.newDoubleRule(WEIGHT),
            AttributeRule.newBooleanRule(RANDOM_SCAN, true),
            new ElementRule(WORKING_PRIOR, new XMLSyntaxRule[]{
                    new ElementRule(FullyConjugateMultivariateTraitLikelihood.class)
            }, true),
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

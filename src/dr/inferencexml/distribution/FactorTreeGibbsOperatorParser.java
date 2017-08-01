package dr.inferencexml.distribution;

import dr.evomodel.continuous.FullyConjugateMultivariateTraitLikelihood;
import dr.evomodel.continuous.GibbsSampleFromTreeInterface;
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
    public final String PATH_PARAMETER = "pathParameter";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        String weightTemp = (String) xo.getAttribute(WEIGHT);
        double weight = Double.parseDouble(weightTemp);
        LatentFactorModel lfm= (LatentFactorModel) xo.getChild(LatentFactorModel.class);
        GibbsSampleFromTreeInterface tree = (GibbsSampleFromTreeInterface) xo.getChild(GibbsSampleFromTreeInterface.class);
        GibbsSampleFromTreeInterface workingTree = null;
        if(xo.getChild(WORKING_PRIOR) != null){
            workingTree = (GibbsSampleFromTreeInterface) xo.getChild(WORKING_PRIOR).getChild(GibbsSampleFromTreeInterface.class);
        }
        boolean randomScan = xo.getAttribute(RANDOM_SCAN, true);

        FactorTreeGibbsOperator lfmOp = new FactorTreeGibbsOperator(weight, lfm, tree, randomScan);
        if(xo.hasAttribute(PATH_PARAMETER)){
            System.out.println("WARNING: Setting Path Parameter is intended for debugging purposes only!");
            lfmOp.setPathParameter(xo.getDoubleAttribute(PATH_PARAMETER));
        }

        return lfmOp;
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(LatentFactorModel.class),
            new ElementRule(GibbsSampleFromTreeInterface.class),
            AttributeRule.newDoubleRule(WEIGHT),
            AttributeRule.newBooleanRule(RANDOM_SCAN, true),
            new ElementRule(WORKING_PRIOR, new XMLSyntaxRule[]{
                    new ElementRule(GibbsSampleFromTreeInterface.class)
            }, true),
            AttributeRule.newDoubleRule(PATH_PARAMETER, true),
    };


    @Override
    public String getParserDescription() {
        return "Gibbs sample a factor row (tip) on a tree";
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

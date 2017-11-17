package dr.inferencexml.operators;

import dr.inference.distribution.IndependentInverseGammaDistributionModel;
import dr.inference.distribution.ShrinkageGibbsOperator;
import dr.inference.model.Parameter;
import dr.xml.*;

public class ShrinkageGibbsOperatorParser extends AbstractXMLObjectParser{
    private final static String SHRINKAGE_GIBBS_OPERATOR = "shrinkageGibbsOperator";
    private final static String LOCAL_PRIOR = "localPrior";
    private final static String GLOBAL_PRIOR = "globalPrior";
    private final static String DATA = "data";
    private final static String WEIGHT = "weight";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        IndependentInverseGammaDistributionModel localPrior = (IndependentInverseGammaDistributionModel) xo.getChild(LOCAL_PRIOR).getChild(IndependentInverseGammaDistributionModel.class);
        IndependentInverseGammaDistributionModel globalPrior = (IndependentInverseGammaDistributionModel) xo.getChild(GLOBAL_PRIOR).getChild(IndependentInverseGammaDistributionModel.class);
        Parameter data = (Parameter) xo.getChild(DATA).getChild(Parameter.class);

        double weight = xo.getDoubleAttribute(WEIGHT);


        return new ShrinkageGibbsOperator(weight, localPrior, globalPrior, data);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            AttributeRule.newDoubleRule(WEIGHT),
            new ElementRule(LOCAL_PRIOR, new XMLSyntaxRule[]{
                    new ElementRule(IndependentInverseGammaDistributionModel.class)
            }),
            new ElementRule(GLOBAL_PRIOR, new XMLSyntaxRule[]{
                    new ElementRule(IndependentInverseGammaDistributionModel.class)
            }),
            new ElementRule(DATA, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            })
    };

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return ShrinkageGibbsOperator.class;
    }

    @Override
    public String getParserName() {
        return SHRINKAGE_GIBBS_OPERATOR;
    }
}

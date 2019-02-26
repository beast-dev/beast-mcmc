package dr.inferencexml.operators.shrinkage;

import dr.inference.distribution.IndependentInverseGammaDistributionModel;
import dr.inference.model.Parameter;
import dr.inference.operators.shrinkage.BayesianBridgeShrinkageOperator;
import dr.xml.*;

public class BayesianBridgeShrinkageOperatorParser extends AbstractXMLObjectParser {

    private final static String BAYESIAN_BRIDGE_PARSER = "shrinkageGibbsOperator";
    private final static String LOCAL_PRIOR = "localPrior";
    private final static String GLOBAL_PRIOR = "globalPrior";
    private final static String DATA = "data";
    private final static String WEIGHT = "weight";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        return new BayesianBridgeShrinkageOperator();
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
        return BayesianBridgeShrinkageOperator.class;
    }

    @Override
    public String getParserName() {
        return BAYESIAN_BRIDGE_PARSER;
    }
}

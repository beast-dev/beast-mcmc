package dr.evomodelxml.operators;

import dr.evomodel.clock.RateEvolutionLikelihood;
import dr.evomodel.operators.RateSampleOperator;
import dr.evomodel.tree.TreeModel;
import dr.inference.operators.MCMCOperator;
import dr.xml.*;

/**
 */
public class RateSampleOperatorParser extends AbstractXMLObjectParser {

    public static final String SAMPLE_OPERATOR = "rateSampleOperator";
    public static final String SAMPLE_ALL = "sampleAll";

    public String getParserName() {
        return SAMPLE_OPERATOR;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {


        final double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);

        final boolean sampleAll = xo.getBooleanAttribute(SAMPLE_ALL);

        TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);

        RateEvolutionLikelihood rateEvolution = (RateEvolutionLikelihood) xo.getChild(RateEvolutionLikelihood.class);

        RateSampleOperator operator = new RateSampleOperator(treeModel, sampleAll, rateEvolution);
        operator.setWeight(weight);
        return operator;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element returns a rateSample operator on a given parameter.";
    }

    public Class getReturnType() {
        return RateSampleOperator.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
            AttributeRule.newBooleanRule(SAMPLE_ALL, true),
            new ElementRule(TreeModel.class),
            new ElementRule(RateEvolutionLikelihood.class, true),
    };
}

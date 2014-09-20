package dr.inferencexml.operators;

import dr.inference.model.Likelihood;
import dr.inference.operators.DirtyLikelihoodOperator;
import dr.inference.operators.MCMCOperator;
import dr.xml.*;

/**
 *
 */
public class DirtyLikelihoodOperatorParser extends AbstractXMLObjectParser {

    public static final String TOUCH_OPERATOR = "dirtyLikelihood";

    public String getParserName() {
        return TOUCH_OPERATOR;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);

        Likelihood likelihood = (Likelihood) xo.getChild(Likelihood.class);

        return new DirtyLikelihoodOperator(likelihood, weight);

    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element returns a operator that forces the entire model likelihood recomputation";
    }

    public Class getReturnType() {
        return DirtyLikelihoodOperator.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
            new ElementRule(Likelihood.class),
    };
}

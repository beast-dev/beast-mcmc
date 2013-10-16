package dr.evomodelxml.clock;

import dr.evomodel.clock.RateEvolutionLikelihood;
import dr.evomodel.clock.UCLikelihood;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 */
public class UCLikelihoodParser extends AbstractXMLObjectParser {

    public static final String UC_LIKELIHOOD = "UCLikelihood";

    public static final String VARIANCE = "variance";

    public String getParserName() {
        return UC_LIKELIHOOD;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);

        Parameter ratesParameter = (Parameter) xo.getElementFirstChild(RateEvolutionLikelihood.RATES);

        Parameter rootRate = (Parameter) xo.getElementFirstChild(RateEvolutionLikelihood.ROOTRATE);

        Parameter variance = (Parameter) xo.getElementFirstChild(VARIANCE);


        boolean isLogSpace = xo.getAttribute(RateEvolutionLikelihood.LOGSPACE, false);

        return new UCLikelihood(tree, ratesParameter, variance, rootRate, isLogSpace);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return
                "This element returns an object that can calculate the likelihood " +
                        "of rates in a tree under the assumption of " +
                        "(log)normally distributed rates. ";
    }

    public Class getReturnType() {
        return UCLikelihood.class;
    }


    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(TreeModel.class),
            new ElementRule(RateEvolutionLikelihood.RATES, Parameter.class, "The branch rates parameter", false),
            AttributeRule.newBooleanRule(RateEvolutionLikelihood.LOGSPACE, true, "true if model considers the log of the rates."),
            new ElementRule(RateEvolutionLikelihood.ROOTRATE, Parameter.class, "The root rate parameter"),
            new ElementRule(VARIANCE, Parameter.class, "The standard deviation of the (log)normal distribution"),
    };

}

package dr.evomodelxml.clock;

import dr.evomodel.clock.ACLikelihood;
import dr.evomodel.clock.RateEvolutionLikelihood;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 */
public class ACLikelihoodParser extends AbstractXMLObjectParser {

    public static final String AC_LIKELIHOOD = "ACLikelihood";

    public static final String VARIANCE = "variance";
    public static final String SHAPE = "shape";

    public static final String DISTRIBUTION = "distribution";

    public String getParserName() {
        return AC_LIKELIHOOD;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);

        Parameter ratesParameter = (Parameter) xo.getElementFirstChild(RateEvolutionLikelihood.RATES);

        Parameter rootRate = (Parameter) xo.getElementFirstChild(RateEvolutionLikelihood.ROOTRATE);

        Parameter variance = (Parameter) xo.getElementFirstChild(VARIANCE);

        boolean isEpisodic = xo.getBooleanAttribute(RateEvolutionLikelihood.EPISODIC);

        //Distribution distributionModel = new InverseGaussianDistribution(0,1);
        //Parameter distribution = (Parameter) xo.getElementFirstChild(DISTRIBUTION);
        String distribution = xo.getStringAttribute(DISTRIBUTION);

        //boolean isLogSpace = xo.getAttribute(LOGSPACE, false);

        //return new ACLikelihood(tree, ratesParameter, variance, rootRate, isEpisodic, isLogSpace);
        return new ACLikelihood(tree, ratesParameter, variance, rootRate, isEpisodic, distribution);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return
                "This element returns an object that can calculate the likelihood " +
                        "of rate changes in a tree under the assumption of " +
                        "distributed rate changes among lineages. " +
                        //"(log)normally distributed rate changes among lineages. " +
                        "Specifically, each branch is assumed to draw a rate from a " +
                        "distribution with mean of the rate in the " +
                        //"(log)normal distribution with mean of the rate in the " +
                        "parent branch and the given standard deviation (the variance can be optionally proportional to " +
                        "branch length).";
    }

    public Class getReturnType() {
        return ACLikelihood.class;
    }


    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(TreeModel.class),
            new ElementRule(RateEvolutionLikelihood.RATES, Parameter.class, "The branch rates parameter", false),
            AttributeRule.newBooleanRule(RateEvolutionLikelihood.EPISODIC, false, "true if model is branch length independent, false if length-dependent."),
            AttributeRule.newStringRule(DISTRIBUTION, false, "The distribution to use"),
            //AttributeRule.newBooleanRule(LOGSPACE, true, "true if model considers the log of the rates."),
            new ElementRule(RateEvolutionLikelihood.ROOTRATE, Parameter.class, "The root rate parameter"),
            new ElementRule(VARIANCE, Parameter.class, "The standard deviation of the distribution"),
            //new ElementRule(VARIANCE, Parameter.class, "The standard deviation of the (log)normal distribution"),
            //new ElementRule(DISTRIBUTION, Parameter.class, "The distribution to use"),
            //new ElementRule(DISTRIBUTION, ParametricDistributionModel.class, "The distribution model for rates among branches", false),
    };
}

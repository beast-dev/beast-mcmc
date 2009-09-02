package dr.evomodel.clock;

import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.math.distributions.NormalDistribution;
import dr.math.MathUtils;
import dr.xml.*;

/**
 * Calculates the likelihood of a set of rate changes in a tree, assuming a (log)normally distributed
 * change in rate at each node, with a mean of the previous (log) rate and a variance proportional to branch length.
 * cf Yang and Rannala 2006
 *
 * @author Michael Defoin Platel
 */
public class ACLikelihood extends RateEvolutionLikelihood {

    public static final String AC_LIKELIHOOD = "ACLikelihood";

    public static final String VARIANCE = "variance";

    public ACLikelihood(TreeModel tree, Parameter ratesParameter, Parameter variance, Parameter rootRate, boolean isEpisodic, boolean isLogSpace) {

        super((isLogSpace) ? "LogNormally Distributed" : "Normally Distributed", tree, ratesParameter, rootRate, isEpisodic);

        this.isLogSpace = isLogSpace;
        this.variance = variance;

        addVariable(variance);
    }

    /**
     * @return the log likelihood of the rate change from the parent to the child.
     */
    double branchRateChangeLogLikelihood(double parentRate, double childRate, double time) {
        double var = variance.getParameterValue(0);

        if (!isEpisodic())
            var *= time;

        if (isLogSpace) {
            double logParentRate = Math.log(parentRate);
            double logChildRate = Math.log(childRate);

            return NormalDistribution.logPdf(logChildRate, logParentRate - (var / 2.), Math.sqrt(var)) - logChildRate;

        } else {
            return NormalDistribution.logPdf(childRate, parentRate, Math.sqrt(var));
        }
    }

    double branchRateSample(double parentRate, double time) {

        //System.out.println(parentRate);

        double var = variance.getParameterValue(0);

        if (!isEpisodic())
            var *= time;

        if (isLogSpace) {
            final double logParentRate = Math.log(parentRate);

            return Math.exp(MathUtils.nextGaussian() * Math.sqrt(var) + logParentRate - (var / 2.));
        } else {
            return MathUtils.nextGaussian() * Math.sqrt(var) + parentRate;
        }
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return AC_LIKELIHOOD;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);

            Parameter ratesParameter = (Parameter) xo.getElementFirstChild(RATES);

            Parameter rootRate = (Parameter) xo.getElementFirstChild(ROOTRATE);

            Parameter variance = (Parameter) xo.getElementFirstChild(VARIANCE);

            boolean isEpisodic = xo.getBooleanAttribute(EPISODIC);

            boolean isLogSpace = xo.getAttribute(LOGSPACE, false);

            return new ACLikelihood(tree, ratesParameter, variance, rootRate, isEpisodic, isLogSpace);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return
                    "This element returns an object that can calculate the likelihood " +
                            "of rate changes in a tree under the assumption of " +
                            "(log)normally distributed rate changes among lineages. " +
                            "Specifically, each branch is assumed to draw a rate from a " +
                            "(log)normal distribution with mean of the rate in the " +
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
                new ElementRule(RATES, Parameter.class, "The branch rates parameter", false),
                AttributeRule.newBooleanRule(EPISODIC, false, "true if model is branch length independent, false if length-dependent."),
                AttributeRule.newBooleanRule(LOGSPACE, true, "true if model considers the log of the rates."),
                new ElementRule(ROOTRATE, Parameter.class, "The root rate parameter"),
                new ElementRule(VARIANCE, Parameter.class, "The standard deviation of the (log)normal distribution"),
        };
    };

    private Parameter variance;
    boolean isLogSpace = false;
}

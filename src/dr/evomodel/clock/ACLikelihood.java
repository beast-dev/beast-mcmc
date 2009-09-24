package dr.evomodel.clock;

import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.math.distributions.NormalDistribution;
import dr.math.distributions.InverseGaussianDistribution;
import dr.math.distributions.LogNormalDistribution;
import dr.math.MathUtils;
import dr.xml.*;

/**
 * Calculates the likelihood of a set of rate changes in a tree, assuming a (log)normal or inverse gaussian distributed
 * change in rate at each node, with a mean of the previous (log) rate and a variance proportional to branch length.
 * cf Yang and Rannala 2006
 *
 * @author Michael Defoin Platel
 * @author Wai Lok Sibon Li
 */
public class ACLikelihood extends RateEvolutionLikelihood {

    public static final String AC_LIKELIHOOD = "ACLikelihood";

    public static final String VARIANCE = "variance";
    public static final String SHAPE = "shape";

    public static final String DISTRIBUTION = "distribution";
    public static final String LOGNORMAL = "logNormal";
    public static final String NORMAL = "normal";
    public static final String INVERSEGAUSSIAN = "inverseGaussian";

    public ACLikelihood(TreeModel tree, Parameter ratesParameter, Parameter variance, Parameter rootRate, boolean isEpisodic, String distribution) {

        //super((isLogSpace) ? "LogNormally Distributed" : "Normally Distributed", tree, ratesParameter, rootRate, isEpisodic);
        super(distribution, tree, ratesParameter, rootRate, isEpisodic);

        //this.isLogSpace = isLogSpace;
        this.variance = variance;
        this.distribution = distribution;

        addVariable(variance);
    }

    /**
     * @return the log likelihood of the rate change from the parent to the child.
     */
    double branchRateChangeLogLikelihood(double parentRate, double childRate, double time) {
        double var = variance.getParameterValue(0);

        if (!isEpisodic())
            var *= time;

        //if (isLogSpace) {
        //    double logParentRate = Math.log(parentRate);
        //    double logChildRate = Math.log(childRate);

        //    return NormalDistribution.logPdf(logChildRate, logParentRate - (var / 2.), Math.sqrt(var)) - logChildRate;

        //} else {
        //    return NormalDistribution.logPdf(childRate, parentRate, Math.sqrt(var));
        //}
        if(distribution.equals(LOGNORMAL)) {
            return LogNormalDistribution.logPdf(childRate, Math.log(parentRate) - (var / 2.), Math.sqrt(var));
        }
        else if(distribution.equals(NORMAL)) {
            return NormalDistribution.logPdf(childRate, parentRate, Math.sqrt(var));
        }
        else if(distribution.equals(INVERSEGAUSSIAN)) { /* Inverse Gaussian */
            double shape = (parentRate * parentRate * parentRate) / var;
            return InverseGaussianDistribution.logPdf(childRate, parentRate, shape);
        }
        else {
            throw new RuntimeException ("Parameter for distribution is not recognised");
        }
    }

    double branchRateSample(double parentRate, double time) {

        double var = variance.getParameterValue(0);

        if (!isEpisodic())
            var *= time;

        //if (isLogSpace) {
        //    final double logParentRate = Math.log(parentRate);

        //    return Math.exp(MathUtils.nextGaussian() * Math.sqrt(var) + logParentRate - (var / 2.));
        //} else {
        //    return MathUtils.nextGaussian() * Math.sqrt(var) + parentRate;
        //}

        if(distribution.equals(LOGNORMAL)) {
            final double logParentRate = Math.log(parentRate);

            return Math.exp(MathUtils.nextGaussian() * Math.sqrt(var) + logParentRate - (var / 2.));
        }
        else if(distribution.equals(NORMAL)) {
            return MathUtils.nextGaussian() * Math.sqrt(var) + parentRate;
        }
        else { /* Inverse Gaussian */
            //return Math.random()
            //Random rand = new Random();
            double lambda = (parentRate * parentRate * parentRate) / var;
            return MathUtils.nextInverseGaussian(parentRate, lambda);
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
                new ElementRule(RATES, Parameter.class, "The branch rates parameter", false),
                AttributeRule.newBooleanRule(EPISODIC, false, "true if model is branch length independent, false if length-dependent."),
                AttributeRule.newStringRule(DISTRIBUTION, false, "The distribution to use"), 
                //AttributeRule.newBooleanRule(LOGSPACE, true, "true if model considers the log of the rates."),
                new ElementRule(ROOTRATE, Parameter.class, "The root rate parameter"),
                new ElementRule(VARIANCE, Parameter.class, "The standard deviation of the distribution"),
                //new ElementRule(VARIANCE, Parameter.class, "The standard deviation of the (log)normal distribution"),
                //new ElementRule(DISTRIBUTION, Parameter.class, "The distribution to use"),
                //new ElementRule(DISTRIBUTION, ParametricDistributionModel.class, "The distribution model for rates among branches", false),
        };
    };

    private Parameter variance;
    //boolean isLogSpace = false;
    String distribution;
}

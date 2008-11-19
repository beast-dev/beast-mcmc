package dr.evomodel.clock;

import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.math.distributions.NormalDistribution;
import dr.math.MathUtils;
import dr.xml.*;

/**
 * Calculates the likelihood of a set of rate changes in a tree, assuming that rates are lognormally distributed
 * cf Yang and Rannala 2006
 *
 * @author Michael Defoin Platel
 */
public class UCLikelihood extends RateEvolutionLikelihood {

    public static final String UC_LIKELIHOOD = "UCLikelihood";
    public static final String ROOTRATE = "rootRate";
    public static final String VARIANCE = "variance";

    public UCLikelihood(TreeModel tree, Parameter ratesParameter, Parameter variance, Parameter rootRate, boolean isLogSpace) {

        super((isLogSpace) ? "LogNormally Distributed" : "Normally Distributed", tree, ratesParameter, rootRate, false);

        this.isLogSpace = isLogSpace;
        this.variance = variance;

        addParameter(variance);
    }

    /**
     * @return the log likelihood of the rate.
     */
    double branchRateChangeLogLikelihood(double foo1, double rate, double foo2) {
        double var = variance.getParameterValue(0);
        double meanRate = rootRateParameter.getParameterValue(0);


        if (isLogSpace) {
            final double logmeanRate = Math.log(meanRate);
            final double logRate = Math.log(rate);

            return NormalDistribution.logPdf(logRate, logmeanRate - (var / 2.), Math.sqrt(var)) - logRate;

        } else {
            return NormalDistribution.logPdf(rate, meanRate, Math.sqrt(var));
        }
    }

    double branchRateSample(double foo1, double foo2) {
        double meanRate = rootRateParameter.getParameterValue(0);

        double var = variance.getParameterValue(0);

        if (isLogSpace) {
            final double logMeanRate = Math.log(meanRate);

            return Math.exp(MathUtils.nextGaussian() * Math.sqrt(var) + logMeanRate - (var / 2.));
        } else {
            return MathUtils.nextGaussian() * Math.sqrt(var) + meanRate;
        }
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return UC_LIKELIHOOD;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);

            Parameter ratesParameter = (Parameter) xo.getElementFirstChild(RATES);

            Parameter rootRate = (Parameter) xo.getElementFirstChild(ROOTRATE);

            Parameter variance = (Parameter) xo.getElementFirstChild(VARIANCE);


            boolean isLogSpace = xo.getAttribute(LOGSPACE, false);

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

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new ElementRule(TreeModel.class),
                new ElementRule(RATES, Parameter.class, "The branch rates parameter", false),
                AttributeRule.newBooleanRule(LOGSPACE, true, "true if model considers the log of the rates."),
                new ElementRule(ROOTRATE, Parameter.class, "The root rate parameter"),
                new ElementRule(VARIANCE, Parameter.class, "The standard deviation of the (log)normal distribution"),
        };
    };

    private Parameter variance;

    boolean isLogSpace = false;
}

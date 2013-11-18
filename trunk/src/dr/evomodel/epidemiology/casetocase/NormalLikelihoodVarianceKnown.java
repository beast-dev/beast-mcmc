package dr.evomodel.epidemiology.casetocase;

import dr.inference.distribution.NormalDistributionModel;
import dr.inference.model.Parameter;
import dr.math.distributions.Distribution;
import dr.math.distributions.NormalDistribution;
import dr.xml.*;

/**
 * Created with IntelliJ IDEA.
 * User: mhall
 * Date: 05/11/2013
 * Time: 15:02
 * To change this template use File | Settings | File Templates.
 */
public class NormalLikelihoodVarianceKnown extends AnalyticallySolvablePosteriorFunction {

    public static final String NORMAL_LIKELIHOOD_VARIANCE_KNOWN = "normalLikelihoodVarianceKnown";

    private Parameter meanPriorMean;
    private Parameter meanPriorVariance;
    private Parameter variance;

    private final NormalDistributionModel meanPrior;
    private final NormalDistributionModel priorPredictive;

    public NormalLikelihoodVarianceKnown(Parameter variance, Parameter meanPriorMean, Parameter meanPriorVariance){
        this.meanPriorMean = meanPriorMean;
        this.meanPriorVariance = meanPriorVariance;
        this.variance = variance;


        meanPrior = new NormalDistributionModel(meanPriorMean, meanPriorVariance, false);
        priorPredictive = new NormalDistributionModel(meanPriorMean,
                new Parameter.Default(Math.sqrt(variance.getParameterValue(0)
                        + meanPriorVariance.getParameterValue(0))));
    }

    // prior probability of the mean taking this value

    public double getPriorValue(Parameter params){
        return meanPrior.pdf(params.getParameterValue(0));
    }

    // probability that a data point will take this value using prior information only

    public double getPriorPredictivePDF(double value){
        return priorPredictive.pdf(value);
    }

    public double getPriorPredictiveCDF(double value){
        return priorPredictive.cdf(value);
    }

    public double getPriorPredictiveInterval(double start, double end){
        return getPriorPredictiveCDF(end) - getPriorPredictiveCDF(start);
    }

    public Distribution getPriorPredictiveDistribution(){
        return priorPredictive;
    }

    public double getPosteriorValue(Parameter params, double[] data){
        double mpv = meanPriorVariance.getParameterValue(0);
        double var = variance.getParameterValue(0);

        int n = data.length;

        double sum = meanPriorMean.getParameterValue(0)/mpv;

        for (double aData : data) {
            sum += aData / var;
        }

        double postVar = 1/((1/mpv) + (n/var));

        double postMean = sum*postVar;

        return NormalDistribution.pdf(params.getParameterValue(0), postMean, postVar);
    }

    public double getLogPosteriorValue(Parameter params, double[] data){
        return Math.log(getPosteriorValue(params, data));
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public static final String MEAN_PRIOR_MEAN = "meanPriorMean";
        public static final String MEAN_PRIOR_VARIANCE = "meanPriorVar";
        public static final String KNOWN_DATA_VARIANCE = "knownDataVar";

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            Parameter mpm = (Parameter)xo.getElementFirstChild(MEAN_PRIOR_MEAN);
            Parameter mpv = (Parameter)xo.getElementFirstChild(MEAN_PRIOR_VARIANCE);
            Parameter kdv = (Parameter)xo.getElementFirstChild(KNOWN_DATA_VARIANCE);

            return new NormalLikelihoodVarianceKnown(kdv, mpm, mpv);
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        public String getParserDescription() {
            return "This element represents a posterior distribution function where the likelihood is normally " +
                    "distributed with known variance and the prior on the mean is the appropriate " +
                    "normally-distributed conjugate";
        }

        public Class getReturnType() {
            return NormalLikelihoodVarianceKnown.class;
        }

        public String getParserName() {
            return NORMAL_LIKELIHOOD_VARIANCE_KNOWN;
        }

        private final XMLSyntaxRule[] rules = {
                new ElementRule(MEAN_PRIOR_MEAN, Parameter.class,
                        "The mean of the normal prior on the mean of the data values"),
                new ElementRule(MEAN_PRIOR_VARIANCE, Parameter.class,
                        "The variance of the normal prior on the mean of the data values"),
                new ElementRule(KNOWN_DATA_VARIANCE, Parameter.class, "The known variance of the data values")
        };

    };

}

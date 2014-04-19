package dr.evomodel.epidemiology.casetocase.periodpriors;

import dr.inference.loggers.LogColumn;
import dr.inference.model.Parameter;
import dr.math.distributions.NormalDistribution;
import dr.math.functionEval.GammaFunction;
import dr.xml.*;

import java.util.ArrayList;
import java.util.Arrays;

/**
 The assumption here is that the periods are drawn from a normal distribution with unknown mean and variance.
 The hyperprior is the conjugate, normal-gamma distribution.

 @author Matthew Hall
 */
public class KnownVarianceNormalPeriodPriorDistribution extends AbstractPeriodPriorDistribution {

    public static final String NORMAL = "knownVarianceNormalPeriodPriorDistribution";
    public static final String LOG = "log";
    public static final String ID = "id";

    // This gets confusing. The data is assumed to be normally distributed with mean mu and stdev sigma. Sigma is
    // known. The prior on mu is that it is _also_ normally distributed with mean mu_0 and stdev sigma_0.

    public static final String MU_0 = "mu0";
    public static final String SIGMA = "sigma";
    public static final String SIGMA_0 = "sigma0";

    private NormalDistribution hyperprior;

    private Parameter posteriorMean;
    private Parameter posteriorVariance;
    private double sigma;

    public KnownVarianceNormalPeriodPriorDistribution(String name, boolean log, double sigma,
                                                      NormalDistribution hyperprior){
        super(name, log);
        this.hyperprior = hyperprior;
        posteriorVariance = new Parameter.Default(1);
        posteriorMean = new Parameter.Default(1);
        addVariable(posteriorVariance);
        addVariable(posteriorMean);
        this.sigma = sigma;
    }

    public KnownVarianceNormalPeriodPriorDistribution(String name, boolean log, double sigma,
                                                      double mu_0, double sigma_0){
        this(name, log, sigma, new NormalDistribution(mu_0, sigma_0));
    }

    public double calculateLogLikelihood(double[] values){

        int count = values.length;

        double mu_0 = hyperprior.getMean();
        double sigma_0 = hyperprior.getSD();

        double var = Math.pow(sigma, 2);
        double var_0 = Math.pow(sigma_0, 2);

        double sum = 0;
        double sumOfSquares = 0;
        for (Double infPeriod : values) {
            sum += infPeriod;
            sumOfSquares += Math.pow(infPeriod, 2);
        }
        double mean = sum/count;

        posteriorMean.setParameterValue(0, ((mu_0/var_0) + sum/var)/(1/var_0 + count/var));
        posteriorVariance.setParameterValue(0, 1/(1/var_0 + count/var));


        logL = Math.log(sigma)
                - count * Math.log(Math.sqrt(2*Math.PI)*sigma)
                - Math.log(Math.sqrt(count*var_0 + var))
                + -sumOfSquares/(2*var) - Math.pow(mu_0, 2)/(2*var_0)
                + (Math.pow(sigma_0*count*mean/sigma, 2) + Math.pow(sigma*mu_0/sigma_0, 2) + 2*count*mean*mu_0)
                /(2*(count*var_0 + var));


        return logL;
    }

    public LogColumn[] getColumns() {
        ArrayList<LogColumn> columns = new ArrayList<LogColumn>(Arrays.asList(super.getColumns()));

        columns.add(new LogColumn.Abstract(getModelName()+"_posteriorMean"){
            protected String getFormattedValue() {
                return String.valueOf(posteriorMean.getParameterValue(0));
            }
        });

        columns.add(new LogColumn.Abstract(getModelName()+"_posteriorVariance"){
            protected String getFormattedValue() {
                return String.valueOf(posteriorVariance.getParameterValue(0));
            }
        });

        return columns.toArray(new LogColumn[columns.size()]);
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return NORMAL;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            String id = (String) xo.getAttribute(ID);

            double mu_0 = xo.getDoubleAttribute(MU_0);
            double sigma = xo.getDoubleAttribute(SIGMA);
            double sigma_0 = xo.getDoubleAttribute(SIGMA_0);

            boolean log;
            log = xo.hasAttribute(LOG) ? xo.getBooleanAttribute(LOG) : false;

            return new KnownVarianceNormalPeriodPriorDistribution(id, log, sigma, mu_0, sigma_0);

        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newBooleanRule(LOG, true),
                AttributeRule.newStringRule(ID, false),
                AttributeRule.newDoubleRule(MU_0, false),
                AttributeRule.newDoubleRule(SIGMA, false),
                AttributeRule.newDoubleRule(SIGMA_0, false),
        };

        public String getParserDescription() {
            return "Calculates the probability of a set of doubles being drawn from the prior posterior distribution" +
                    "of a normal distribution of unknown mean and known standard deviation sigma";
        }

        public Class getReturnType() {
            return KnownVarianceNormalPeriodPriorDistribution.class;
        }
    };

}

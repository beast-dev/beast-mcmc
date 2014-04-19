package dr.evomodel.epidemiology.casetocase.periodpriors;

import dr.evomodel.epidemiology.casetocase.periodpriors.AbstractPeriodPriorDistribution;
import dr.inference.loggers.LogColumn;
import dr.inference.model.Parameter;
import dr.math.distributions.NormalGammaDistribution;
import dr.math.functionEval.GammaFunction;
import dr.xml.*;

import java.util.ArrayList;
import java.util.Arrays;

/**
 The assumption here is that the periods are drawn from a normal distribution with unknown mean and variance.
 The hyperprior is the conjugate, normal-gamma distribution.

 @author Matthew Hall
 */
public class NormalPeriodPriorDistribution extends AbstractPeriodPriorDistribution {

    public static final String NORMAL = "normalPeriodPriorDistribution";
    public static final String LOG = "log";
    public static final String ID = "id";
    public static final String MU = "mu";
    public static final String LAMBDA = "lambda";
    public static final String ALPHA = "alpha";
    public static final String BETA = "beta";

    private NormalGammaDistribution hyperprior;

    private Parameter posteriorMean;
    private Parameter posteriorBeta;
    private Parameter posteriorExpectedPrecision;

    public NormalPeriodPriorDistribution(String name, boolean log, NormalGammaDistribution hyperprior){
        super(name, log);
        this.hyperprior = hyperprior;
        posteriorBeta = new Parameter.Default(1);
        posteriorMean = new Parameter.Default(1);
        posteriorExpectedPrecision = new Parameter.Default(1);
        addVariable(posteriorBeta);
        addVariable(posteriorMean);
        addVariable(posteriorExpectedPrecision);
    }

    public NormalPeriodPriorDistribution(String name, boolean log, double mu_0, double lambda_0,
                                         double alpha_0, double beta_0){
        this(name, log, new NormalGammaDistribution(mu_0, lambda_0, alpha_0, beta_0));
    }



    public double calculateLogLikelihood(double[] values){

        int count = values.length;

        double[] infPredictiveDistributionParameters=hyperprior.getParameters();

        double mu_0 = infPredictiveDistributionParameters[0];
        double lambda_0 = infPredictiveDistributionParameters[1];
        double alpha_0 = infPredictiveDistributionParameters[2];
        double beta_0 = infPredictiveDistributionParameters[3];

        double lambda_n = lambda_0 + count;
        double alpha_n = alpha_0 + count/2;
        double sum = 0;
        for (Double infPeriod : values) {
            sum += infPeriod;
        }
        double mean = sum/count;

        double sumOfDifferences = 0;
        for (Double infPeriod : values) {
            sumOfDifferences += Math.pow(infPeriod-mean,2);
        }

        posteriorMean.setParameterValue(0, (lambda_0*mu_0 + sum)/(lambda_0 + count));

        double beta_n = beta_0 + 0.5*sumOfDifferences
                + lambda_0*count*Math.pow(mean-mu_0, 2)/(2*(lambda_0+count));

        posteriorBeta.setParameterValue(0, beta_n);
        posteriorExpectedPrecision.setParameterValue(0, alpha_n/beta_n);

        logL = GammaFunction.logGamma(alpha_n)
                - GammaFunction.logGamma(alpha_0)
                + alpha_0*Math.log(beta_0)
                - alpha_n*Math.log(beta_n)
                + 0.5*Math.log(lambda_0)
                - 0.5*Math.log(lambda_n)
                - (count/2)*Math.log(2*Math.PI);


        return logL;
    }

    public LogColumn[] getColumns() {
        ArrayList<LogColumn> columns = new ArrayList<LogColumn>(Arrays.asList(super.getColumns()));

        columns.add(new LogColumn.Abstract(getModelName()+"_posteriorMean"){
            protected String getFormattedValue() {
                return String.valueOf(posteriorMean.getParameterValue(0));
            }
        });

        columns.add(new LogColumn.Abstract(getModelName()+"_posteriorBeta"){
            protected String getFormattedValue() {
                return String.valueOf(posteriorBeta.getParameterValue(0));
            }
        });

        columns.add(new LogColumn.Abstract(getModelName()+"_posteriorExpectedPrecision"){
            protected String getFormattedValue() {
                return String.valueOf(posteriorExpectedPrecision.getParameterValue(0));
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

            double mu = xo.getDoubleAttribute(MU);
            double lambda = xo.getDoubleAttribute(LAMBDA);
            double alpha = xo.getDoubleAttribute(ALPHA);
            double beta = xo.getDoubleAttribute(BETA);

            boolean log;
            log = xo.hasAttribute(LOG) ? xo.getBooleanAttribute(LOG) : false;

            return new NormalPeriodPriorDistribution(id, log, mu, lambda, alpha, beta);

        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newBooleanRule(LOG, true),
                AttributeRule.newStringRule(ID, false),
                AttributeRule.newDoubleRule(MU, false),
                AttributeRule.newDoubleRule(LAMBDA, false),
                AttributeRule.newDoubleRule(ALPHA, false),
                AttributeRule.newDoubleRule(BETA, false)
        };

        public String getParserDescription() {
            return "Calculates the probability of a set of doubles being drawn from the prior posterior distribution" +
                    "of a normal distribution of unknown mean and variance";
        }

        public Class getReturnType() {
            return NormalPeriodPriorDistribution.class;
        }
    };

}

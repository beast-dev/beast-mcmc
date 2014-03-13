package dr.evomodel.epidemiology.casetocase.periodpriors;

import dr.evomodel.epidemiology.casetocase.periodpriors.AbstractPeriodPriorDistribution;
import dr.math.distributions.ExponentialDistribution;
import dr.xml.*;
import org.apache.commons.math.stat.descriptive.DescriptiveStatistics;

/**
 * This isn't terribly rigorous, truth be told
 *
 * Created by mhall on 28/02/2014.
 */
public class ExponentiallyDistributedStDevPeriodPriorDistribution extends AbstractPeriodPriorDistribution {

    public static final String EXPONENTIALLY_DISTRIBUTED_STDEV = "exponentiallyDistributedStDevPeriodPriorDistribution";
    public static final String LOG = "log";
    public static final String ID = "id";
    public static final String LAMBDA = "lambda";

    private ExponentialDistribution hyperprior;


    public ExponentiallyDistributedStDevPeriodPriorDistribution(String name, boolean log,
                                                                ExponentialDistribution hyperprior){
        super(name, log);
        this.hyperprior = hyperprior;
    }

    public ExponentiallyDistributedStDevPeriodPriorDistribution(String name, boolean log,
                                                                double lambda){
        this(name, log, new ExponentialDistribution(lambda));
    }

    public double calculateLogLikelihood(double[] values){

        DescriptiveStatistics stats = new DescriptiveStatistics(values);

        logL = hyperprior.logPdf(stats.getStandardDeviation());

        return logL;

    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return EXPONENTIALLY_DISTRIBUTED_STDEV;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            String id = (String) xo.getAttribute(ID);

            boolean log;
            log = xo.hasAttribute(LOG) ? xo.getBooleanAttribute(LOG) : false;

            double lambda = xo.getDoubleAttribute(LAMBDA);

            return new ExponentiallyDistributedStDevPeriodPriorDistribution(id, log, lambda);

        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newBooleanRule(LOG, true),
                AttributeRule.newStringRule(ID, false),
                AttributeRule.newDoubleRule(LAMBDA, false)
        };

        public String getParserDescription() {
            return "Calculates the probability of observing a list of doubles under the assumption that their " +
                    "standard deviation is exponentially distributed";
        }

        public Class getReturnType() {
            return ExponentiallyDistributedStDevPeriodPriorDistribution.class;
        }
    };

}

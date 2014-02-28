package dr.evomodel.epidemiology.casetocase.PeriodPriors;

import dr.xml.*;
import org.apache.commons.math.stat.descriptive.DescriptiveStatistics;

/**
 * Created by mhall on 28/02/2014.
 */
public class OneOverStDevInfectiousOrLatentPeriodPriorDistribution extends InfectiousOrLatentPeriodPriorDistribution {

    public static final String ONE_OVER_STDEV = "oneOverStDevInfectiousOrLatentPeriodPriorDistribution";
    public static final String LOG = "log";
    public static final String ID = "id";

    public OneOverStDevInfectiousOrLatentPeriodPriorDistribution(String name, boolean log){
        super(name, log);
    }

    public double calculateLogLikelihood(double[] values){

        DescriptiveStatistics stats = new DescriptiveStatistics(values);

        logL.setParameterValue(0,-Math.log(stats.getStandardDeviation()));

        return -Math.log(stats.getStandardDeviation());

    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return ONE_OVER_STDEV;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            String id = (String) xo.getAttribute(ID);

            boolean log;
            log = xo.hasAttribute(LOG) ? xo.getBooleanAttribute(LOG) : false;

            return new OneOverStDevInfectiousOrLatentPeriodPriorDistribution(id, log);

        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newBooleanRule(LOG, true),
                AttributeRule.newStringRule(ID, false)
        };

        public String getParserDescription() {
            return "Calculates the probability of observing a list of doubles with probability proportional to" +
                    "1 over their standard deviation (the Jeffreys prior for normally distributed data)";
        }

        public Class getReturnType() {
            return OneOverStDevInfectiousOrLatentPeriodPriorDistribution.class;
        }
    };

}

package dr.evomodel.epidemiology.casetocase.periodpriors;

import dr.xml.*;
import org.apache.commons.math.stat.descriptive.DescriptiveStatistics;

/**
 * Probability is proportion to one over variance
 */
public class OneOverVariancePeriodPriorDistribution extends AbstractPeriodPriorDistribution {

    public static final String ONE_OVER_VARIANCE = "oneOverVariancePeriodPriorDistribution";
    public static final String LOG = "log";
    public static final String ID = "id";

    public OneOverVariancePeriodPriorDistribution(String name, boolean log){
        super(name, log);
    }

    public double calculateLogLikelihood(double[] values){

        DescriptiveStatistics stats = new DescriptiveStatistics(values);

        logL = -Math.log(stats.getVariance());

        return logL;

    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return ONE_OVER_VARIANCE;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            String id = (String) xo.getAttribute(ID);

            boolean log;
            log = xo.hasAttribute(LOG) ? xo.getBooleanAttribute(LOG) : false;

            return new OneOverVariancePeriodPriorDistribution(id, log);

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
                    "1 over their variance";
        }

        public Class getReturnType() {
            return OneOverVariancePeriodPriorDistribution.class;
        }
    };



}

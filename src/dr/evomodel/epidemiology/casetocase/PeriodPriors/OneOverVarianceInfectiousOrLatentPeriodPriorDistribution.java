package dr.evomodel.epidemiology.casetocase.PeriodPriors;

import dr.xml.*;
import org.apache.commons.math.stat.descriptive.DescriptiveStatistics;

/**
 * Created by mhall on 28/02/2014.
 */
public class OneOverVarianceInfectiousOrLatentPeriodPriorDistribution extends InfectiousOrLatentPeriodPriorDistribution {

    public static final String ONE_OVER_VARIANCE = "oneOverVarianceInfectiousOrLatentPeriodPriorDistribution";
    public static final String LOG = "log";
    public static final String ID = "id";

    public OneOverVarianceInfectiousOrLatentPeriodPriorDistribution(String name, boolean log){
        super(name, log);
    }

    public double calculateLogLikelihood(double[] values){

        DescriptiveStatistics stats = new DescriptiveStatistics(values);

        logL.setParameterValue(0,-Math.log(stats.getStandardDeviation()));

        return -Math.log(stats.getVariance());

    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return ONE_OVER_VARIANCE;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            String id = (String) xo.getAttribute(ID);

            boolean log;
            log = xo.hasAttribute(LOG) ? xo.getBooleanAttribute(LOG) : false;

            return new OneOverVarianceInfectiousOrLatentPeriodPriorDistribution(id, log);

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
            return OneOverVarianceInfectiousOrLatentPeriodPriorDistribution.class;
        }
    };



}

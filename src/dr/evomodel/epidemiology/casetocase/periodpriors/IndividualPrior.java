package dr.evomodel.epidemiology.casetocase.periodpriors;

import dr.evomodel.epidemiology.casetocase.AbstractCase;
import dr.inference.distribution.ParametricDistributionModel;
import dr.xml.*;

import java.util.HashMap;

/**
 * Created by mhall on 03/06/2014.
 */
public class IndividualPrior extends AbstractPeriodPriorDistribution {

    ParametricDistributionModel distribution;

    public static final String INDIVIDUAL_PRIOR = "individualPrior";
    public static final String ID = "id";
    public static final String DISTRIBUTION = "distribution";

    public IndividualPrior(String name, boolean log, ParametricDistributionModel distribution){
        super(name, log);
        this.distribution = distribution;
    }

    public double calculateLogLikelihood(double[] values) {
        double out = 0;

        for(int i=0; i<values.length; i++){
            out += distribution.pdf(values[i]);
        }

        return out;
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return INDIVIDUAL_PRIOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            String id = (String) xo.getAttribute(ID);

            return null;

        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newStringRule(ID, false),
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

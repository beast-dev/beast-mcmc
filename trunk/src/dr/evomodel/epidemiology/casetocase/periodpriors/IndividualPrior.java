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

    public IndividualPrior(String name, ParametricDistributionModel distribution){
        super(name, false);
        this.distribution = distribution;
    }

    public double calculateLogLikelihood(double[] values) {
        double out = 0;

        for (double value : values) {
            out += distribution.logPdf(value);
        }

        return out;
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return INDIVIDUAL_PRIOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            String id = (String) xo.getAttribute(ID);
            ParametricDistributionModel distribution =
                    (ParametricDistributionModel)xo.getElementFirstChild(DISTRIBUTION);

            return new IndividualPrior(id, distribution);

        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newStringRule(ID, false),
                new ElementRule(DISTRIBUTION, ParametricDistributionModel.class)
        };

        public String getParserDescription() {
            return "Calculates the probability of a set of doubles all being drawn from the specified prior " +
                    "distribution";
        }

        public Class getReturnType() {
            return IndividualPrior.class;
        }
    };
}

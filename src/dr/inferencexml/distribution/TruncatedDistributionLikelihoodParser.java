package dr.inferencexml.distribution;

import dr.inference.distribution.DistributionLikelihood;
import dr.inference.model.Parameter;
import dr.inference.model.TruncatedDistributionLikelihood;
import dr.xml.*;

/**
 * @author Max Tolkoff
 */
public class TruncatedDistributionLikelihoodParser extends AbstractXMLObjectParser {
    public static final String TRUNCATED_DISTRIBUTION_LIKELIHOOD = "truncatedDistributionLikelihood";
    public static final String LOW = "low";
    public static final String HIGH = "high";

    @Override
    public String getParserName() {
        return TRUNCATED_DISTRIBUTION_LIKELIHOOD;
    }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        DistributionLikelihood likelihood = (DistributionLikelihood) xo.getChild(DistributionLikelihood.class);
        Parameter low;
        Parameter high;
        if(xo.getChild(LOW) != null){
            low = (Parameter) xo.getChild(LOW).getChild(Parameter.class);
        }
        else{
            low = new Parameter.Default("low", Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY);
        }
        if(xo.getChild(HIGH) != null){
            high = (Parameter) xo.getChild(HIGH).getChild(Parameter.class);
        }
        else{
            high = new Parameter.Default("high", Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY);
        }

        return new TruncatedDistributionLikelihood(likelihood, low, high);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    @Override
    public String getParserDescription() {
        return "Produces a truncated distribution likelihood";
    }

    @Override
    public Class getReturnType() {
        return TruncatedDistributionLikelihood.class;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(DistributionLikelihood.class),
            new ElementRule(LOW, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }, true),
            new ElementRule(HIGH, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }, true),
    };
}

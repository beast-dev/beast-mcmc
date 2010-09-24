package dr.inferencexml.distribution;

import dr.inference.distribution.MixedDistributionLikelihood;
import dr.inference.distribution.ParametricDistributionModel;
import dr.inference.model.Statistic;
import dr.xml.*;

/**
 * Reads a distribution likelihood from a DOM Document element.
 */
public class MixedDistributionLikelihoodParser extends AbstractXMLObjectParser {
    public static final String DISTRIBUTION_LIKELIHOOD = "mixedDistributionLikelihood";

    public static final String DISTRIBUTION0 = "distribution0";
    public static final String DISTRIBUTION1 = "distribution1";
    public static final String DATA = "data";
    public static final String INDICATORS = "indicators";

    public String getParserName() {
        return DISTRIBUTION_LIKELIHOOD;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        XMLObject cxo0 = xo.getChild(DISTRIBUTION0);
        ParametricDistributionModel model0 = (ParametricDistributionModel) cxo0.getChild(ParametricDistributionModel.class);

        XMLObject cxo1 = xo.getChild(DISTRIBUTION1);
        ParametricDistributionModel model1 = (ParametricDistributionModel) cxo1.getChild(ParametricDistributionModel.class);

        Statistic data = (Statistic) ((XMLObject) xo.getChild(DATA)).getChild(Statistic.class);
        Statistic indicators = (Statistic) ((XMLObject) xo.getChild(INDICATORS)).getChild(Statistic.class);

        ParametricDistributionModel[] models = {model0, model1};
        try {
          return new MixedDistributionLikelihood(models, data, indicators);
        } catch( Exception e) {
            throw new XMLParseException(e.getMessage());
        }
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(DISTRIBUTION0,
                    new XMLSyntaxRule[]{new ElementRule(ParametricDistributionModel.class)}),
            new ElementRule(DISTRIBUTION1,
                    new XMLSyntaxRule[]{new ElementRule(ParametricDistributionModel.class)}),
            new ElementRule(DATA, new XMLSyntaxRule[]{new ElementRule(Statistic.class)}),
            new ElementRule(INDICATORS, new XMLSyntaxRule[]{new ElementRule(Statistic.class)}),
    };

    public String getParserDescription() {
        return "Calculates the likelihood of some data given some mix of parametric distributions.";
    }

    public Class getReturnType() {
        return MixedDistributionLikelihood.class;
    }

}

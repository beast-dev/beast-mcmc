package dr.evomodelxml.coalescent;

import dr.evomodel.coalescent.MultilocusNPCoalescentLikelihoodGradient;
import dr.evomodel.coalescent.MultilocusNonparametricCoalescentLikelihood;
import dr.inference.model.Parameter;
import dr.xml.*;

public class MultilocusNPCoalescentLikelihoodGradientParser extends AbstractXMLObjectParser {
    private static final String PARSER_NAME = "multilocusNPCoalescentLikelihoodGradient";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        MultilocusNonparametricCoalescentLikelihood likelihood = (MultilocusNonparametricCoalescentLikelihood) xo.getChild(MultilocusNonparametricCoalescentLikelihood.class);
        Parameter logPopSizes = (Parameter) xo.getChild(Parameter.class);
        return new MultilocusNPCoalescentLikelihoodGradient(likelihood, logPopSizes);
    }

    public XMLSyntaxRule[] getSyntaxRules() {return rules;}

    private final XMLSyntaxRule[] rules = {
            new ElementRule(MultilocusNonparametricCoalescentLikelihood.class),
            new ElementRule(Parameter.class),
    };

    @Override
    public String getParserDescription() {
        return "Gradient of the likelihood of a multilocus nonparametric coalescent model" +
                "with respect to the log population sizes";
    }

    @Override
    public Class getReturnType() {
        return MultilocusNPCoalescentLikelihoodGradient.class;
    }

    @Override
    public String getParserName() {
        return PARSER_NAME;
    }
}

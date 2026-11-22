package dr.evomodelxml.coalescent.nonparametric;

import dr.evomodel.coalescent.nonparametric.NonParametricCoalescentLikelihoodGradient;
import dr.evomodel.coalescent.nonparametric.NonParametricCoalescentLikelihood;
import dr.inference.model.Parameter;
import dr.xml.*;

public class NonParametricCoalescentLikelihoodGradientParser extends AbstractXMLObjectParser {
    private static final String PARSER_NAME = "NonParametricCoalescentLikelihoodGradient";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        NonParametricCoalescentLikelihood likelihood = (NonParametricCoalescentLikelihood) xo.getChild(NonParametricCoalescentLikelihood.class);
        Parameter logPopSizes = (Parameter) xo.getChild(Parameter.class);
        return new NonParametricCoalescentLikelihoodGradient(likelihood, logPopSizes);
    }

    public XMLSyntaxRule[] getSyntaxRules() {return rules;}

    private final XMLSyntaxRule[] rules = {
            new ElementRule(NonParametricCoalescentLikelihood.class),
            new ElementRule(Parameter.class),
    };

    @Override
    public String getParserDescription() {
        return "Gradient of the likelihood of a multilocus nonparametric coalescent model" +
                "with respect to the log population sizes";
    }

    @Override
    public Class getReturnType() {
        return NonParametricCoalescentLikelihoodGradient.class;
    }

    @Override
    public String getParserName() {
        return PARSER_NAME;
    }
}

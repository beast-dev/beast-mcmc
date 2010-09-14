package dr.inferencexml.model;

import dr.inference.model.Likelihood;
import dr.inference.model.ThreadedCompoundLikelihood;
import dr.xml.*;

/**
 *
 */
public class ThreadedCompoundLikelihoodParser extends AbstractXMLObjectParser {

    public static final String THREADED_COMPOUND_LIKELIHOOD = "threadedCompoundLikelihood";
    public static final String WEIGHT = "robustWeight";

    public String getParserName() {
        return THREADED_COMPOUND_LIKELIHOOD;
    }

    public String[] getParserNames() {
        return new String[]{getParserName()};
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        ThreadedCompoundLikelihood compoundLikelihood = new ThreadedCompoundLikelihood();

        for (int i = 0; i < xo.getChildCount(); i++) {
            if (xo.getChild(i) instanceof Likelihood) {
                compoundLikelihood.addLikelihood((Likelihood) xo.getChild(i));
            } else {

                Object rogueElement = xo.getChild(i);

                throw new XMLParseException("An element (" + rogueElement + ") which is not a likelihood has been added to a " + THREADED_COMPOUND_LIKELIHOOD + " element");
            }
        }

        double weight = xo.getAttribute(WEIGHT, 0.0);
        if (weight < 0)
            throw new XMLParseException("Robust weight must be non-negative.");
        compoundLikelihood.setWeightFactor(Math.exp(-weight));

        return compoundLikelihood;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A likelihood function which is simply the product of its component likelihood functions.";
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(Likelihood.class, 1, Integer.MAX_VALUE),
            AttributeRule.newDoubleRule(WEIGHT, true),
    };

    public Class getReturnType() {
        return ThreadedCompoundLikelihood.class;
    }
}

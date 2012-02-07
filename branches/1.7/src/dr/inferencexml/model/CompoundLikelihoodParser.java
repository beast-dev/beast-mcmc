package dr.inferencexml.model;

import dr.inference.model.CompoundLikelihood;
import dr.inference.model.Likelihood;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

/**
 */
public class CompoundLikelihoodParser extends AbstractXMLObjectParser {
    public static final String COMPOUND_LIKELIHOOD = "compoundLikelihood";
    public static final String THREADS = "threads";
    public static final String POSTERIOR = "posterior";
    public static final String PRIOR = "prior";
    public static final String LIKELIHOOD = "likelihood";

    public String getParserName() {
        return COMPOUND_LIKELIHOOD;
    }

    public String[] getParserNames() {
        return new String[]{getParserName(), POSTERIOR, PRIOR, LIKELIHOOD};
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        // the default is 0 threads but an XML attribute can override it
        int threads = xo.getAttribute(THREADS, 0);

//        if (xo.getName().equalsIgnoreCase(LIKELIHOOD)) {
            // if this is '<likelihood>' then the default is to use a cached thread pool...
            threads = xo.getAttribute(THREADS, -1);

            // both the XML attribute and a system property can override it
            if (System.getProperty("thread.count") != null) {

                threads = Integer.parseInt(System.getProperty("thread.count"));
                if (threads < -1 || threads > 1000) {
                    // put an upper limit here - may be unnecessary?
                    threads = -1;
                }
            }
//        }

        List<Likelihood> likelihoods = new ArrayList<Likelihood>();
        for (int i = 0; i < xo.getChildCount(); i++) {
            final Object child = xo.getChild(i);
            if (child instanceof Likelihood) {
                likelihoods.add((Likelihood) child);
            } else {

                throw new XMLParseException("An element (" + child + ") which is not a likelihood has been added to a "
                        + COMPOUND_LIKELIHOOD + " element");
            }
        }

        CompoundLikelihood compoundLikelihood;

        if (xo.getName().equalsIgnoreCase(POSTERIOR)) {
            compoundLikelihood = new CompoundLikelihood(threads, likelihoods);
        } else {
            compoundLikelihood = new CompoundLikelihood(0, likelihoods);
        }

        if (compoundLikelihood.getThreadCount() > 0) {
            Logger.getLogger("dr.evomodel").info("Likelihood is using " + threads + " threads.");
//            } else if (threads < 0) {
//                Logger.getLogger("dr.evomodel").info("Likelihood is using a cached thread pool.");
        }

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
            AttributeRule.newIntegerRule(THREADS, true),
            new ElementRule(Likelihood.class, 1, Integer.MAX_VALUE)
    };

    public Class getReturnType() {
        return CompoundLikelihood.class;
    }
}

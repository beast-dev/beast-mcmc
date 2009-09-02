package dr.app.beagle.evomodel.parsers;

import dr.inference.model.*;
import dr.xml.*;

import java.util.logging.Logger;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @author Marc Suchard
 */
public class CompoundLikelihoodParser extends AbstractXMLObjectParser {

    public static final String LIKELIHOOD = "likelihood";


    public String getParserName() {
        return LIKELIHOOD;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
   
        ThreadedCompoundLikelihood compoundLikelihood = new ThreadedCompoundLikelihood();

        for (int i = 0; i < xo.getChildCount(); i++) {
            if (xo.getChild(i) instanceof Likelihood) {
                compoundLikelihood.addLikelihood((Likelihood) xo.getChild(i));
            } else {

                Object rogueElement = xo.getChild(i);

                throw new XMLParseException("An element (" + rogueElement + ") which is not a likelihood has been added to the " + LIKELIHOOD + " element");
            }
        }

        Logger.getLogger("dr.evomodel").info("Multithreaded Likelihood, using " + compoundLikelihood.getLikelihoodCount() + " threads.");

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
    };

    public Class getReturnType() {
        return CompoundLikelihood.class;
    }
}
package dr.evomodel.beagle.parsers;

import dr.inference.model.Parameter;
import dr.inference.model.Likelihood;
import dr.inference.model.CompoundModel;
import dr.inference.model.CompoundLikelihood;
import dr.xml.*;
import dr.evomodel.beagle.substmodel.FrequencyModel;
import dr.evomodel.beagle.substmodel.HKY;

import java.util.logging.Logger;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Callable;
import java.util.ArrayList;
import java.util.List;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 */
public class CompoundLikelihoodParser extends AbstractXMLObjectParser {

    public static final String LIKELIHOOD = "likelihood";


    public String getParserName() {
        return LIKELIHOOD;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        int threads = 0;
        for (int i = 0; i < xo.getChildCount(); i++) {
            if (xo.getChild(i) instanceof Likelihood) {
                threads ++;
            }
        }

        CompoundLikelihood compoundLikelihood = new CompoundLikelihood(threads);

        for (int i = 0; i < xo.getChildCount(); i++) {
            if (xo.getChild(i) instanceof Likelihood) {
                compoundLikelihood.addLikelihood((Likelihood) xo.getChild(i));
            } else {

                Object rogueElement = xo.getChild(i);

                throw new XMLParseException("An element (" + rogueElement + ") which is not a likelihood has been added to the " + LIKELIHOOD + " element");
            }
        }

        Logger.getLogger("dr.evomodel").info("Multithreaded Likelihood, using " + threads + " threads.");


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
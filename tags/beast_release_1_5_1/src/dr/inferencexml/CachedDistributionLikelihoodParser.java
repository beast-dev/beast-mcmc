package dr.inferencexml;

import dr.inference.distribution.DistributionLikelihood;
import dr.inference.model.Variable;
import dr.xml.*;

import java.util.logging.Logger;

/**
 * @author Marc Suchard
 */
public class CachedDistributionLikelihoodParser extends AbstractXMLObjectParser {

    public static final String MODEL_NAME = "cachedPrior";

   // public static final String RATE_BLOCK = "rates";
    //public static final String INDICATOR_BLOCK = "indicators";

    public String getParserName() {
        return MODEL_NAME;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        final String name = xo.hasId() ? xo.getId() : MODEL_NAME;

        final DistributionLikelihood likelihood = (DistributionLikelihood) xo.getChild(DistributionLikelihood.class);
        final Variable variable = (Variable) xo.getChild(Variable.class);

        final Logger logger = Logger.getLogger("dr.inference");
        logger.info("Constructing a cache around likelihood '" + likelihood.getId() + "', signal = " + variable.getVariableName());

        return new dr.inference.distribution.CachedDistributionLikelihood(name, likelihood, variable);

    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(DistributionLikelihood.class),
            new ElementRule(Variable.class),
    };

    public String getParserDescription() {
        return "Calculates a cached likelihood of some data given some parametric or empirical distribution.";
    }

    public Class getReturnType() {
        return DistributionLikelihood.class;
    }
}

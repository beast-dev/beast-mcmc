package dr.inferencexml.distribution;

import dr.inference.distribution.DistributionLikelihood;
import dr.inference.distribution.TwoPartsDistributionLikelihood;
import dr.inference.model.Parameter;
import dr.xml.*;
import dr.math.distributions.Distribution;

/**
 * @author Chieh-Hsi Wu
 *
 * Parser for TwoPartDistribution likelihood
 */
public class TwoPartsDistributionLikelihoodParser extends AbstractXMLObjectParser{
    public static final String TWO_PART_DISTRIBUTION_LIKELIHOOD = "twoPartDistribution";
    public static final String PRIOR = "priorLik";
    public static final String PSEUDO_PRIOR = "pseudoPriorLik";
    public static final String PARAMETER_VECTOR = "parameterVector";
    public static final String PARAMETER_INDEX = "paramIndex";
    public static final String SELECTED_VARIABLE = "selectedVariable";


    public String getParserName() {
        return TWO_PART_DISTRIBUTION_LIKELIHOOD;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        DistributionLikelihood priorLikelihood = (DistributionLikelihood)xo.getElementFirstChild(PRIOR);
        DistributionLikelihood pseudoPriorLikelihood = (DistributionLikelihood)xo.getElementFirstChild(PSEUDO_PRIOR);
        Distribution prior = priorLikelihood.getDistribution();
        Distribution pseudoPrior = pseudoPriorLikelihood.getDistribution();
        Parameter bitVector = (Parameter)xo.getElementFirstChild(PARAMETER_VECTOR);
        int paramIndex = xo.getIntegerAttribute(PARAMETER_INDEX);
        Parameter selectedVariable = (Parameter)xo.getElementFirstChild(SELECTED_VARIABLE);


        TwoPartsDistributionLikelihood likelihood =
                new TwoPartsDistributionLikelihood(
                        prior,
                        pseudoPrior,
                        bitVector,
                        paramIndex
                );
        likelihood.addData(selectedVariable);

        return likelihood;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {

    };

    public String getParserDescription() {
        return "Calculates the likelihood of some data given some parametric or empirical distribution.";
    }

    public Class getReturnType() {
        return TwoPartsDistributionLikelihood.class;
    }
}

package dr.evomodelxml.epidemiology;

import dr.evomodel.epidemiology.TwoPathogenModel;
import dr.evomodel.epidemiology.TwoPathogenModelBirthRates;
import dr.evomodel.epidemiology.TwoPathogenModelLikelihood;
import dr.inference.model.Parameter;
import dr.xml.*;

public class TwoPathogenModelLikelihoodParser extends AbstractXMLObjectParser {

    public static final String TWO_PATHOGEN_MODEL_LIKELIHOOD = "TwoPathogenModelLikelihood";
    //public static final String TWO_PATHOGEN_MODEL = "TwoPathogenModel";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        TwoPathogenModel tpm = (TwoPathogenModel) xo.getChild(TwoPathogenModel.class);

        return new TwoPathogenModelLikelihood(tpm);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element returns stochastic trajectory likelihood for a two-pathogen compartmental model (Shrestha et al., 2011)" +
                "that is used to parameterize birth-death process priors. ";
    }

    public Class getReturnType() {
        return TwoPathogenModelLikelihood.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    public String getParserName() {
        return TWO_PATHOGEN_MODEL_LIKELIHOOD;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(TwoPathogenModel.class),
    };

}
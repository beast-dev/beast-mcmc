package dr.inferencexml.model;

import dr.inference.model.IntermediateReturnLikelihood;
import dr.inference.model.IntermediateReturnParameter;
import dr.inference.model.LatentFactorModel;
import dr.xml.*;

/**
 * Created by max on 10/22/14.
 */
public class IntermediateReturnParameterParser extends AbstractXMLObjectParser {
    final static String INTERMEDIATE_RETURN_PARAMETER = "IntermediateReturnParameter";
    private final XMLSyntaxRule[] rules = {
            new ElementRule(LatentFactorModel.class),
    };

    @Override
    public Object parseXMLObject (XMLObject xo)throws XMLParseException {
        return new IntermediateReturnParameter((IntermediateReturnLikelihood) xo.getChild(LatentFactorModel.class)) {
        };
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules () {
        return rules;
    }

    @Override
    public String getParserDescription () {
        return "Gets Latent Factor Model to return data with residuals computed";
    }

    @Override
    public Class getReturnType () {
        return IntermediateReturnLikelihood.class;
    }

    @Override
    public String getParserName () {
        return INTERMEDIATE_RETURN_PARAMETER;
    }
};

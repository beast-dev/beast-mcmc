package dr.inferencexml.model;

import dr.inference.model.LFMParameter;
import dr.inference.model.LatentFactorModel;
import dr.xml.*;

/**
 * Created by max on 10/22/14.
 */
public class LFMParameterParser extends AbstractXMLObjectParser {
    final static String LFM_PARAMETER = "LFMParameter";
    private final XMLSyntaxRule[] rules = {
            new ElementRule(LatentFactorModel.class),
    };

    @Override
    public Object parseXMLObject (XMLObject xo)throws XMLParseException {
        return new LFMParameter((LatentFactorModel) xo.getChild(LatentFactorModel.class));
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
        return LFMParameter.class;
    }

    @Override
    public String getParserName () {
        return LFM_PARAMETER;
    }
};

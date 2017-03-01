package dr.inferencexml.model;

import dr.inference.model.LFMFactorPotentialDerivative;
import dr.inference.model.LatentFactorModel;
import dr.xml.*;

/**
 * @author Max Tolkoff
 */
public class LogLFMFactorPotentialDerivativeParser extends AbstractXMLObjectParser {
    public final static String LFM_FACTOR_DERIVATIVE = "logLFMFactorDerivative";

    @Override
    public String getParserName() {
        return LFM_FACTOR_DERIVATIVE;
    }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        LatentFactorModel lfm = (LatentFactorModel) xo.getChild(LatentFactorModel.class);

        return new LFMFactorPotentialDerivative(lfm);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(LatentFactorModel.class),
    };

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return LFMFactorPotentialDerivative.class;
    }
}

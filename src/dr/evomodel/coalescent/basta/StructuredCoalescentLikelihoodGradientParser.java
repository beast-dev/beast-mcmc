package dr.evomodel.coalescent.basta;

import dr.evomodel.substmodel.SubstitutionModel;
import dr.xml.*;

public class StructuredCoalescentLikelihoodGradientParser extends AbstractXMLObjectParser {

    private static final String NAME = "structuredCoalescentLikelihoodGradient";
    private static final String WRT_PARAMETER = "wrtParameter";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        BastaLikelihood BastaLikelihood = (BastaLikelihood) xo.getChild(BastaLikelihood.class);
        SubstitutionModel substitutionModel = (SubstitutionModel) xo.getChild(SubstitutionModel.class);
        String wrtParamter = (String) xo.getAttribute(WRT_PARAMETER);

        StructuredCoalescentLikelihoodGradient.WrtParameter type = StructuredCoalescentLikelihoodGradient.WrtParameter.factory(wrtParamter);

        return new StructuredCoalescentLikelihoodGradient(BastaLikelihood, substitutionModel, type);
//        return StructuredCoalescentLikelihoodGradient.factory(BastaLikelihood, substitutionModel, type);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newStringRule(WRT_PARAMETER),
            new ElementRule(BastaLikelihood.class),
            new ElementRule(SubstitutionModel.class),
    };

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return StructuredCoalescentLikelihoodGradient.class;
    }

    @Override
    public String getParserName() {
        return NAME;
    }
}

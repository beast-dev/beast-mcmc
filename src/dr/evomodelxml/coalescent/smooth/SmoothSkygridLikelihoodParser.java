package dr.evomodelxml.coalescent.smooth;

import dr.evomodel.coalescent.smooth.SmoothSkygridLikelihood;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

public class SmoothSkygridLikelihoodParser extends AbstractXMLObjectParser {

    private static final String PARSER_NAME = "smoothSkygridLikelihood";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        return new SmoothSkygridLikelihood(null, null,
                null , null);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[0];
    }

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return SmoothSkygridLikelihood.class;
    }

    @Override
    public String getParserName() {
        return PARSER_NAME;
    }
}

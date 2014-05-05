package dr.inferencexml.operators;

import dr.math.Polynomial;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

/**
 * Created by Max on 5/5/14.
 */
public class FactorGibbsSamplerParser extends AbstractXMLObjectParser {
    private final String FACTOR_GIBBS_SAMPLER="factorGibbsSampler";
    private final String WEIGHT="weight";


    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        String weightTemp= (String) xo.getAttribute(WEIGHT);
        Double weight=Double.parseDouble(weightTemp);
        return null;
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
        return null;
    }

    @Override
    public String getParserName() {
        return FACTOR_GIBBS_SAMPLER;
    }
}

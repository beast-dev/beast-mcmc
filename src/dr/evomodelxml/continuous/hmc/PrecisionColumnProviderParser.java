
package dr.evomodelxml.continuous.hmc;

import dr.inference.hmc.PrecisionColumnProvider;
import dr.inference.model.MatrixParameter;
import dr.inference.model.MatrixParameterInterface;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

/**
 * @author Zhenyu Zhang
 */

public class PrecisionColumnProviderParser extends AbstractXMLObjectParser {
    private static final String PRODUCT_PROVIDER = "precisionColumn";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        MatrixParameter matrix = (MatrixParameter) xo.getChild(MatrixParameterInterface.class);
        return new PrecisionColumnProvider.Generic(matrix);
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
        return PRODUCT_PROVIDER;
    }
}

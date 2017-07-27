package dr.inferencexml.model;


import dr.inference.model.ElementWiseMatrixMultiplicationParameter;
import dr.inference.model.MatrixParameter;
import dr.inference.model.MatrixParameterInterface;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

/**
 * Created by max on 11/30/15.
 */
public class ElementWiseMatrixMultiplicationParser extends AbstractXMLObjectParser {
    public final static String ELEMENT_WISE_MATRIX_MULTIPLICATION_PARAMETER="elementWiseMatrixMultiplicationParameter";
    public final static String NAME="name";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        final String name = xo.hasId() ? xo.getId() : null;
        MatrixParameterInterface[] matList=new MatrixParameterInterface[xo.getChildCount()];
        for (int i = 0; i <xo.getChildCount(); i++) {
            matList[i]=(MatrixParameterInterface) xo.getChild(i);
        }

        return new ElementWiseMatrixMultiplicationParameter(name, matList);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[0];
    }

    @Override
    public String getParserDescription() {
        return "Returns element wise matrix multiplication of a series of matrices";
    }

    @Override
    public Class getReturnType() {
        return ElementWiseMatrixMultiplicationParameter.class;
    }

    @Override
    public String getParserName() {
        return ELEMENT_WISE_MATRIX_MULTIPLICATION_PARAMETER;
    }
}

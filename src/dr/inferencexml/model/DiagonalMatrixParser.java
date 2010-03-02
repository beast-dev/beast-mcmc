package dr.inferencexml.model;

import dr.inference.model.DiagonalMatrix;
import dr.inference.model.MatrixParameter;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 *
 */
public class DiagonalMatrixParser extends AbstractXMLObjectParser {

    public final static String MATRIX_PARAMETER = "diagonalMatrix";

    public String getParserName() {
        return MATRIX_PARAMETER;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {


        Parameter diagonalParameter = (Parameter) xo.getChild(Parameter.class);

        return new DiagonalMatrix(diagonalParameter);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A diagonal matrix parameter constructed from its diagonals.";
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(Parameter.class, 1, 1),
    };

    public Class getReturnType() {
        return MatrixParameter.class;
    }
    
}

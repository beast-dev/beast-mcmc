package dr.inferencexml.model;

import dr.inference.model.CompoundSymmetricMatrix;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 *
 */
public class CompoundSymmetricMatrixParser extends AbstractXMLObjectParser {

    public final static String MATRIX_PARAMETER = "compoundSymmetricMatrix";
    public static final String DIAGONAL = "diagonal";
    public static final String OFF_DIAGONAL = "offDiagonal";
    public static final String AS_CORRELATION = "asCorrelation";

    public String getParserName() {
        return MATRIX_PARAMETER;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        XMLObject cxo = xo.getChild(DIAGONAL);
        Parameter diagonalParameter = (Parameter) cxo.getChild(Parameter.class);

        cxo = xo.getChild(OFF_DIAGONAL);
        Parameter offDiagonalParameter = (Parameter) cxo.getChild(Parameter.class);

        boolean asCorrelation = xo.getAttribute(AS_CORRELATION, false);

        return new CompoundSymmetricMatrix(diagonalParameter, offDiagonalParameter, asCorrelation);
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
            new ElementRule(DIAGONAL,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(OFF_DIAGONAL,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            AttributeRule.newBooleanRule(AS_CORRELATION, true)
    };

    public Class getReturnType() {
        return CompoundSymmetricMatrix.class;
    }
}

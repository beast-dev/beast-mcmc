package dr.inferencexml.model;

import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Parameter;
import dr.xml.*;

public class MatrixColumnParameterParser extends AbstractXMLObjectParser {

    private static final String MATRIX_COLUMN = "matrixColumnParameter";
    private static final String COLUMN = "column";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        MatrixParameterInterface matParam = (MatrixParameterInterface) xo.getChild(MatrixParameterInterface.class);
        int col = xo.getIntegerAttribute(COLUMN);
        if (col > matParam.getColumnDimension()) {
            throw new XMLParseException("Requested column " + col + ", but the matrix only has " +
                    matParam.getColumnDimension() + " columns.");
        }
        return matParam.getParameter(col - 1);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                new ElementRule(MatrixParameterInterface.class),
                AttributeRule.newIntegerRule(COLUMN)
        };
    }

    @Override
    public String getParserDescription() {
        return "Returns parameter associated with single column of matrix parameter.";
    }

    @Override
    public Class getReturnType() {
        return Parameter.class;
    }

    @Override
    public String getParserName() {
        return MATRIX_COLUMN;
    }
}

package dr.inferencexml.model;

import dr.inference.model.BlockDiagonalMatrixParameter;
import dr.inference.model.MatrixParameterInterface;
import dr.xml.*;

public class BlockDiagonalMatrixParameterParser extends AbstractXMLObjectParser {

    private final static String BLOCK_DIAGONAL_MATRIX_PARAMETER = "blockDiagonalMatrixParameter";
    private final static String REPLICATES = "replicates";

    public String getParserName() {
        return BLOCK_DIAGONAL_MATRIX_PARAMETER;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        String name = xo.hasId() ? xo.getId() : null;
        MatrixParameterInterface block = (MatrixParameterInterface) xo.getChild(MatrixParameterInterface.class);
        int replicates = xo.getIntegerAttribute(REPLICATES);

        return new BlockDiagonalMatrixParameter(name, block, replicates);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A matrix parameter constructed from its component parameters.";
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newIntegerRule(REPLICATES),
            new ElementRule(MatrixParameterInterface.class),
    };

    public Class getReturnType() {
        return BlockDiagonalMatrixParameter.class;
    }
}

package dr.inferencexml.model;

import dr.inference.model.BlockUpperTriangularMatrixParameter;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
@author Max Tolkoff
 */
public class BlockUpperTriangularMatrixParameterParser extends AbstractXMLObjectParser {
    private static final String BLOCK_UPPER_TRIANGULAR_MATRIX="blockUpperTriangularMatrixParameter";
    private static final String COLUMN_DIMENSION="columnDimension";
    private static final String TRANSPOSE="transpose";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        final String name = xo.hasId() ? xo.getId() : null;
        final boolean transpose= xo.getAttribute(TRANSPOSE, false);
//        int rowDim=xo.getChildCount();
//        int colDim;
        Parameter temp=null;
//        if(xo.hasAttribute(COLUMN_DIMENSION)) {
//            colDim = xo.getAttribute(COLUMN_DIMENSION, 1);
//        }
//        else
//        {
//            temp=(Parameter) xo.getChild(xo.getChildCount()-1);
//            colDim=temp.getDimension();
//        }

        Parameter[] params=new Parameter[xo.getChildCount()];



        for (int i = 0; i < xo.getChildCount(); i++) {
            temp = (Parameter) xo.getChild(i);
            params[i]=temp;}

        BlockUpperTriangularMatrixParameter ltmp=new BlockUpperTriangularMatrixParameter(name, params);
        if(transpose){
            return ltmp.transposeBlock();
        }
        else {
            return ltmp;  //To change body of implemented methods use File | Settings | File Templates.
        }
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(Parameter.class, 0, Integer.MAX_VALUE),
            AttributeRule.newBooleanRule(TRANSPOSE, true),
            AttributeRule.newIntegerRule(COLUMN_DIMENSION, true),
    };


    @Override
    public String getParserDescription() {
        return "Returns a blockUpperTriangularMatrixParameter which is a compoundParameter which forces the last element to be of full length, the second to last element to be of full length-1, etc.";  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public Class getReturnType() {
        return BlockUpperTriangularMatrixParameter.class;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public String getParserName() {
        return BLOCK_UPPER_TRIANGULAR_MATRIX;  //To change body of implemented methods use File | Settings | File Templates.
    }
}

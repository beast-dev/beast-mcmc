package dr.inferencexml.model;

import com.sun.tools.internal.xjc.reader.xmlschema.bindinfo.BIConversion;
import dr.inference.model.LowerTriangularMatrixParameter;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
@author Max Tolkoff
 */
public class LowerTriangularMatrixParameterParser extends AbstractXMLObjectParser {
    private static final String LOWER_TRIANGULAR_MATRIX="LowerTriangularMatrixParameter";
    private static final String COLUMN_DIMENSION="columnDimension";


    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        final String name = xo.hasId() ? xo.getId() : null;

        int rowDim=xo.getChildCount();
        int colDim;
        Parameter temp=null;
        if(xo.hasAttribute(COLUMN_DIMENSION)) {
            colDim = xo.getAttribute(COLUMN_DIMENSION, 1);
        }
        else
        {
            temp=(Parameter) xo.getChild(xo.getChildCount()-1);
            colDim=temp.getDimension();
        }

        Parameter[] params=new Parameter[rowDim];



        for (int i = 0; i < xo.getChildCount(); i++) {
            temp = (Parameter) xo.getChild(i);
            params[i]=temp;}

        LowerTriangularMatrixParameter ltmp=new LowerTriangularMatrixParameter(name, params, rowDim, colDim);

        return ltmp;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(Parameter.class, 0, Integer.MAX_VALUE),
            AttributeRule.newIntegerRule(COLUMN_DIMENSION, true),
    };


    @Override
    public String getParserDescription() {
        return "Returns a lowerTriangularMatrixParameter which is a compoundParameter which forces the last element to be of full length, the second to last element to be of full length-1, etc.";  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public Class getReturnType() {
        return LowerTriangularMatrixParameter.class;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public String getParserName() {
        return LOWER_TRIANGULAR_MATRIX;  //To change body of implemented methods use File | Settings | File Templates.
    }
}

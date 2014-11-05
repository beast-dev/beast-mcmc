package dr.inferencexml.model;

import dr.inference.model.MatrixMatrixProduct;
import dr.inference.model.MatrixParameter;
import dr.xml.*;

/**
 * Created by max on 10/22/14.
 */
public class MatrixMatrixProductParser extends AbstractXMLObjectParser {
    final static String MATRIX_MATRIX_PRODUCT = "MatrixMatrixProduct";
    final static String LEFT="left";
    final static String RIGHT="right";
    private final XMLSyntaxRule[] rules = {
            new ElementRule(LEFT, MatrixParameter.class),
            new ElementRule(RIGHT, MatrixParameter.class),
    };

    @Override
    public Object parseXMLObject (XMLObject xo)throws XMLParseException {
        MatrixParameter[] temp=new MatrixParameter[2];
        temp[0]=(MatrixParameter) xo.getChild(LEFT).getChild(MatrixParameter.class);
        temp[1]=(MatrixParameter) xo.getChild(RIGHT).getChild(MatrixParameter.class);
        return new MatrixMatrixProduct(temp) {

        };
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules () {
        return rules;
    }

    @Override
    public String getParserDescription () {
        return "Gets Latent Factor Model to return data with residuals computed";
    }

    @Override
    public Class getReturnType () {
        return MatrixMatrixProduct.class;
    }

    @Override
    public String getParserName () {
        return MATRIX_MATRIX_PRODUCT;
    }
};

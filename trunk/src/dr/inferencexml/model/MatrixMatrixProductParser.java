package dr.inferencexml.model;

import dr.inference.model.CompoundParameter;
import dr.inference.model.MatrixMatrixProduct;
import dr.inference.model.MatrixParameter;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * Created by max on 10/22/14.
 */
public class MatrixMatrixProductParser extends AbstractXMLObjectParser {
    final static String MATRIX_MATRIX_PRODUCT = "MatrixMatrixProduct";
    final static String LEFT="left";
    final static String RIGHT="right";
    final static String IN_PLACE="InPlace";
    final static String COLUMN_MASK="ColumnMask";
    private final XMLSyntaxRule[] rules = {
            new ElementRule(LEFT, MatrixParameter.class),
            new ElementRule(RIGHT, CompoundParameter.class),
            new ElementRule(COLUMN_MASK, Parameter.class, "Only some columns need to be multiplied", true),
            new ElementRule(IN_PLACE, MatrixParameter.class, "Matrix values that are returned", true)
    };

    @Override
    public Object parseXMLObject (XMLObject xo)throws XMLParseException {
        MatrixParameter[] temp=new MatrixParameter[3];
        temp[0]=(MatrixParameter) xo.getChild(LEFT).getChild(MatrixParameter.class);
        temp[1]=MatrixParameter.recast(((CompoundParameter) xo.getChild(RIGHT).getChild(CompoundParameter.class)).getVariableName(), (CompoundParameter) xo.getChild(RIGHT).getChild(CompoundParameter.class));
        if(xo.getChild(IN_PLACE)!=null){
            temp[2]=(MatrixParameter) xo.getChild(IN_PLACE).getChild(MatrixParameter.class);
        }
        else{
            int rowDim=temp[0].getRowDimension();
            int colDim=temp[1].getColumnDimension();
            Parameter[] params=new Parameter[colDim];
            for (int i = 0; i <colDim ; i++) {
                params[i]=new Parameter.Default(rowDim);
            }
            temp[2]=new MatrixParameter(null, params);
        }
        Parameter ColumnMask;
        if(xo.getChild(COLUMN_MASK)!=null)
            ColumnMask=(Parameter) xo.getChild(COLUMN_MASK).getChild(MatrixParameter.class);
        else
            ColumnMask=new Parameter.Default(null, temp[1].getColumnDimension(), 1);
        return new MatrixMatrixProduct(temp, ColumnMask);
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

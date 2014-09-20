package dr.evomodelxml.continuous;

import dr.inference.model.MatrixParameter;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.inference.model.Parameter;
import dr.xml.XMLSyntaxRule;

/**
@author Max Tolkoff
 */

//the existing compound symmetric matrix function creates a matrix that is structurally compound symmetric. This builds
//a compound symmetric matrix that can be non compound symmetric when exposed to a transition kernel
public class BuildCompoundSymmetricMatrix extends AbstractXMLObjectParser {
    public final static String BUILD_DIAGONAL_MATRIX="buildCompoundSymmetricMatrix";
    public final static String DIMENSION="dimension";
    public final static String DIAGONAL="diagonal";
    public final static String OFF_DIAGONAL="offDiagonal";


    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        int dim=xo.getAttribute(DIMENSION,1);
        int diag=xo.getAttribute(DIAGONAL,1);
        int offDiag=xo.getAttribute(OFF_DIAGONAL, 0);
        final String name = xo.hasId() ? xo.getId() : null;
        Parameter[] answer = new Parameter[dim];



        for (int i=0;i<dim; i++){
            answer[i] = new Parameter.Default(dim);
            for(int j=0; j<dim; j++)
        {if(i==j) {
            answer[i].setParameterValue(j, diag);
        }
        else
            {answer[i].setParameterValue(j,offDiag);}
        }}
        return new MatrixParameter(name, answer);
    }

    @Override
    public String getParserName() {
        return BUILD_DIAGONAL_MATRIX;
    }

    @Override
    public Class getReturnType() {
        return MatrixParameter.class;
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[0];  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public String getParserDescription() {
        return "Returns a compound symmetric Matrix Parameter that is flexible, i.e. when exposed to a transition kernel can yield a result that is not compound symmetric";  //To change body of implemented methods use File | Settings | File Templates.
    }
}

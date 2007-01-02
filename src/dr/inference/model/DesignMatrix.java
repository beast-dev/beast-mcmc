package dr.inference.model;

import dr.xml.*;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

/**
 * Created by IntelliJ IDEA.
 * User: msuchard
 * Date: Dec 29, 2006
 * Time: 4:03:54 PM
 * To change this template use File | Settings | File Templates.
 */
public class DesignMatrix extends MatrixParameter {
    public static final String DESIGN_MATRIX = "designMatrix";
    public static final String ADD_INTERCEPT = "addIntercept";

    public DesignMatrix(String name) {
        super(name);
    }

    public DesignMatrix(String name, Parameter[] parameters) {
        super(name, parameters);
    }

    // **************************************************************
    // XMLElement IMPLEMENTATION
    // **************************************************************

    public Element createElement(Document d) {
        throw new RuntimeException("Not implemented yet!");
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {


        public String getParserName() {
            return DESIGN_MATRIX;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            DesignMatrix designMatrix = new DesignMatrix(DESIGN_MATRIX);
            boolean addIntercept = false;


            if (xo.hasAttribute(ADD_INTERCEPT)) {
                addIntercept = xo.getBooleanAttribute(ADD_INTERCEPT);
            }

            int dim = 0;

            for (int i = 0; i < xo.getChildCount(); i++) {
                Parameter parameter = (Parameter) xo.getChild(i);
                designMatrix.addParameter(parameter);
                if (i == 0)
                    dim = parameter.getDimension();
                else if (dim != parameter.getDimension())
                    throw new XMLParseException("All parameters must have the same dimension to construct a rectangular design matrix");
            }

            if (addIntercept) {
                Parameter intercept = new Parameter.Default(dim);
                intercept.setId("Intercept");
                designMatrix.addParameter(intercept);
            }

            return designMatrix;
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

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new ElementRule(Parameter.class, 1, Integer.MAX_VALUE),
        };

        public Class getReturnType() {
            return DesignMatrix.class;
        }
    };
}

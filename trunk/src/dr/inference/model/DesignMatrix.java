package dr.inference.model;

import dr.xml.*;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

/**
 * @author Marc Suchard
 */
public class DesignMatrix extends MatrixParameter {
    public static final String DESIGN_MATRIX = "designMatrix";
    public static final String ADD_INTERCEPT = "addIntercept";
    public static final String FORM = "form";
    public static final String DIMENSION = "dimension";

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
            boolean addIntercept = xo.getAttribute(ADD_INTERCEPT, false);


            int dim = 0;

            if (xo.hasAttribute(FORM)) {
                String type = (String) xo.getStringAttribute(FORM);
                if (type.compareTo("J") == 0) {
                dim = xo.getAttribute(DIMENSION,1);
                System.err.println("dim = "+dim);
                for(int i=0; i<dim; i++) {
                    Parameter parameter = new Parameter.Default(dim);
                    designMatrix.addParameter(parameter);
                }
                } else
                    throw new XMLParseException("Unknown designMatrix form.");
            } else {

            for (int i = 0; i < xo.getChildCount(); i++) {
                Parameter parameter = (Parameter) xo.getChild(i);
                designMatrix.addParameter(parameter);
                if (i == 0)
                    dim = parameter.getDimension();
                else if (dim != parameter.getDimension())
                    throw new XMLParseException("All parameters must have the same dimension to construct a rectangular design matrix");
            }
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
                AttributeRule.newBooleanRule(ADD_INTERCEPT, true),
                new ElementRule(Parameter.class, 0, Integer.MAX_VALUE), // TODO or have the following                            
                AttributeRule.newStringRule(FORM,true),     // TODO Should have to include both FORM and DIMENSION at the same time
                AttributeRule.newIntegerRule(DIMENSION,true)
        };

        public Class getReturnType() {
            return DesignMatrix.class;
        }
    };
}

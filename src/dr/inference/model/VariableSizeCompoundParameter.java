/*
 * VariableSizeCompoundParameter.java
 *
 * (c) 2002-2005 BEAST Development Core Team
 *
 * This package may be distributed under the
 * Lesser Gnu Public Licence (LGPL)
 */

package dr.inference.model;

import dr.xml.*;

/**
 * A multidimensional, variable-size parameter constructed from its component parameters.
 *
 * @author Marc Suchard
 * @version $Id: VariableSizeCompoundParameter.java,v 1.1.1.1 2006/10/02 20:47:46 msuchard Exp $
 */
public class VariableSizeCompoundParameter extends ARGCompoundParameter {

    public static final String COMPOUND_PARAMETER = "variableSizeCompoundParameter";

    public VariableSizeCompoundParameter(String name, Parameter[] parameters) {
        super(name, parameters);
//		this.name = name;
        //this.parameters = parameters;
        //dimension = 0;
        //for (int i =0; i < parameters.length; i++) {
        //	dimension += parameters[i].getDimension();
        //	parameters[i].addParameterListener(this);
        //}
        numParameters = parameters.length;
    }


    public VariableSizeCompoundParameter(String name) {
        super(name);
        numParameters = 0;
        //this.name = name;
        //dimension = 0;
    }

    public Parameter getParameter(int i) {
        return parameters[i];
    }

    public void addParameter(Parameter param) {

        if (parameters == null) {
            parameters = new Parameter[]{param};
        } else {
            Parameter[] newParams = new Parameter[parameters.length + 1];
            System.arraycopy(parameters, 0, newParams, 0, parameters.length);
            newParams[parameters.length] = param;
            parameters = newParams;
        }
        dimension += param.getDimension();
        numParameters++;
        param.addParameterListener(this);
    }

    public void removeParameter(Parameter param) {
        //System.err.println(this.getId()+" contains "+parameters.length+" parameters.");
        //System.err.println("Attempting to remove "+param.getId());
        int len = parameters.length;
        int lenM1 = len - 1;
        Parameter[] newParams = new Parameter[lenM1];
        int index = 0;
        for (int i = 0; i < len; i++) {
            if (parameters[i] != param) {
                //System.err.print("i = "+i+" ");
                //System.err.println("copying at index = "+index);
                if (index < lenM1) {
                    newParams[index] = parameters[i];
                    index++;
                } else
                    throw new IllegalArgumentException("Unable to find " + param.getId() + " in " + getId());
            }
            //} else {
            //	System.err.println("Found at i = "+i);
            //}
        }
        parameters = newParams;
        dimension -= param.getDimension();
        numParameters--;
        param.removeParameterListener(this);
    }

    public int getNumParameters() {
        return numParameters;
    }

    //public void setDimension(int dim) { throw new RuntimeException(); }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return COMPOUND_PARAMETER;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            VariableSizeCompoundParameter compoundParameter = new VariableSizeCompoundParameter(COMPOUND_PARAMETER);

            for (int i = 0; i < xo.getChildCount(); i++) {
                compoundParameter.addParameter((Parameter) xo.getChild(i));
            }
            return compoundParameter;
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "A multidimensional, variable-sized parameter constructed from its component parameters.";
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new ElementRule(Parameter.class, 1, Integer.MAX_VALUE),
        };

        public Class getReturnType() {
            return CompoundParameter.class;
        }
    };

    private int numParameters;
}

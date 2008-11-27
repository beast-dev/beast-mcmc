/*
 * ARGCompoundParameter.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.inference.model;

import dr.xml.*;

/**
 * A multidimensional parameter constructed from its component parameters.
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @version $Id: ARGCompoundParameter.java,v 1.13 2005/06/14 10:40:34 rambaut Exp $
 */
public class ARGCompoundParameter extends Parameter.Abstract implements ParameterListener {

    public static final String COMPOUND_PARAMETER = "compoundParameter";

    public ARGCompoundParameter(String name, Parameter[] parameters) {
//		this.name = name;
        this.parameters = parameters;
        dimension = 0;
        for (Parameter parameter : parameters) {
            dimension += parameter.getDimension();
            parameter.addParameterListener(this);
        }

    }

    public ARGCompoundParameter(String name) {
//		this.name = name;
        dimension = 0;
    }

    public void addParameter(Parameter param) {

        if (parameters == null) {
            parameters = new Parameter[]{param};
        } else {
            Parameter[] newParams = new Parameter[parameters.length + 1];
            for (int i = 0; i < parameters.length; i++) {
                newParams[i] = parameters[i];
            }
            newParams[parameters.length] = param;
            parameters = newParams;
        }
        dimension += param.getDimension();
        param.addParameterListener(this);
    }

    public final String getParameterName() {
        return getId();
    }

    public Parameter getParameter(int index) {
        return parameters[index];
    }

    public int getNumberOfParameters() {
        return parameters.length;
    }

    public final String getDimensionName(int dim) {
        int[] index = new int[1];
        Parameter param = findParameter(dim, index);
        return param.getDimensionName(index[0]);

    }

    public int getDimension() {
        return dimension;
    }

    public void setDimension(int dim) {
        throw new RuntimeException();
    }

    public void addBounds(Bounds boundary) {

        if (bounds == null) createBounds();
        bounds.addBounds(boundary);
    }

    public Bounds getBounds() {

        if (bounds == null) createBounds();
        return bounds;
    }

    public void addDimension(int index, double value) {
        throw new RuntimeException("Not implemented!");
    }

    public double removeDimension(int index) {
        throw new RuntimeException("Not implemented!");
    }

    private void createBounds() {
        bounds = new IntersectionBounds(getDimension());
        bounds.addBounds(new CompoundBounds());
    }

    public double getParameterValue(int dim) {
        int[] index = new int[1];
        Parameter param = findParameter(dim, index);

        return param.getParameterValue(index[0]);
    }

    public double[] inspectParametersValues() {
        throw new RuntimeException("Not yet implemented.");
    }

    public void setParameterValue(int dim, double value) {
        int[] index = new int[1];
        Parameter param = findParameter(dim, index);

        param.setParameterValue(index[0], value);
    }

    public void setParameterValueQuietly(int dim, double value) {
        int[] index = new int[1];
        Parameter param = findParameter(dim, index);

        param.setParameterValueQuietly(index[0], value);
    }

    protected final void storeValues() {
        for (Parameter parameter : parameters) {
            parameter.storeParameterValues();
        }
    }

    protected final void restoreValues() {
        for (Parameter parameter : parameters) {
            parameter.restoreParameterValues();
        }
    }

    protected final void acceptValues() {
        for (Parameter parameter : parameters) {
            parameter.acceptParameterValues();
        }
    }

    protected final void adoptValues(Parameter source) {
        // the parameters that make up a compound parameter will have
        // this function called on them individually so we don't need
        // to do anything here.
    }

    private Parameter findParameter(int dim, int[] outIndex) {
        int k = 0;

        for (Parameter parameter : parameters) {
            if (dim < k + parameter.getDimension()) {
                outIndex[0] = dim - k;
                return parameter;
            }
            k += parameter.getDimension();
        }

        throw new IllegalArgumentException("index out of bound in compound parameter");
    }

    public String toString() {
        StringBuffer buffer = new StringBuffer(String.valueOf(getParameterValue(0)));
        final Bounds bounds = getBounds();
        buffer.append("[").append(String.valueOf(bounds.getLowerLimit(0)));
        buffer.append(",").append(String.valueOf(bounds.getUpperLimit(0))).append("]");

        for (int i = 1; i < getDimension(); i++) {
            buffer.append(", ").append(String.valueOf(getParameterValue(i)));
            buffer.append("[").append(String.valueOf(bounds.getLowerLimit(i)));
            buffer.append(",").append(String.valueOf(bounds.getUpperLimit(i))).append("]");
        }
        return buffer.toString();
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return COMPOUND_PARAMETER;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            ARGCompoundParameter compoundParameter = new ARGCompoundParameter(COMPOUND_PARAMETER);

            for (int i = 0; i < xo.getChildCount(); i++) {
                compoundParameter.addParameter((Parameter) xo.getChild(i));
            }
            return compoundParameter;
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "A multidimensional parameter constructed from its component parameters.";
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new ElementRule(Parameter.class, 1, Integer.MAX_VALUE),
        };

        public Class getReturnType() {
            return ARGCompoundParameter.class;
        }
    };

    // ****************************************************************
    // Parameter listener interface
    // ****************************************************************

    public void parameterChangedEvent(Parameter parameter, int index, Parameter.ChangeType type) {

        int dim = 0;
        for (Parameter parameter1 : parameters) {
            if (parameter == parameter1) {
                fireParameterChangedEvent(dim + index, type);
                break;
            }
            dim += parameter1.getDimension();
        }
    }

    // ****************************************************************
    // Private and protected stuff
    // ****************************************************************

    private class CompoundBounds implements Bounds {

        public double getUpperLimit(int dim) {
            int[] index = new int[1];
            Parameter param = findParameter(dim, index);
            return param.getBounds().getUpperLimit(index[0]);
        }

        public double getLowerLimit(int dim) {
            int[] index = new int[1];
            Parameter param = findParameter(dim, index);
            return param.getBounds().getLowerLimit(index[0]);
        }

        public int getBoundsDimension() {
            return getDimension();
        }
    }


    protected Parameter[] parameters = null;
    private IntersectionBounds bounds = null;
    protected int dimension;
//	private String name;
}

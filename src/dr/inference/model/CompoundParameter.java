/*
 * CompoundParameter.java
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
 * @version $Id: CompoundParameter.java,v 1.13 2005/06/14 10:40:34 rambaut Exp $
 */
public class CompoundParameter extends Parameter.Abstract implements ParameterListener {

    public static final String COMPOUND_PARAMETER = "compoundParameter";

    public CompoundParameter(String name, Parameter[] parameters) {
//		this.name = name;
        this.parameters = parameters;
        dimension = 0;
        for (int i = 0; i < parameters.length; i++) {
            dimension += parameters[i].getDimension();
            parameters[i].addParameterListener(this);
        }

    }

    public CompoundParameter(String name) {
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

    private void createBounds() {
        bounds = new IntersectionBounds(getDimension());
        bounds.addBounds(new CompoundBounds());
    }

    public double getParameterValue(int dim) {
        int[] index = new int[1];
        Parameter param = findParameter(dim, index);

        return param.getParameterValue(index[0]);
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
        for (int i = 0; i < parameters.length; i++) {
            parameters[i].storeParameterValues();
        }
    }

    protected final void restoreValues() {
        for (int i = 0; i < parameters.length; i++) {
            parameters[i].restoreParameterValues();
        }
    }

    protected final void acceptValues() {
        for (int i = 0; i < parameters.length; i++) {
            parameters[i].acceptParameterValues();
        }
    }

    protected final void adoptValues(Parameter source) {
        // the parameters that make up a compound parameter will have
        // this function called on them individually so we don't need
        // to do anything here.
    }

    private Parameter findParameter(int dim, int[] outIndex) {
        int k = 0;

        for (int j = 0; j < parameters.length; j++) {
            if (dim < k + parameters[j].getDimension()) {
                outIndex[0] = dim - k;
                return parameters[j];
            }
            k += parameters[j].getDimension();
        }

        throw new IllegalArgumentException("index out of bound in compound parameter");
    }

    public String toString() {
        StringBuffer buffer = new StringBuffer(String.valueOf(getParameterValue(0)));
        buffer.append("[" + String.valueOf(getBounds().getLowerLimit(0)));
        buffer.append("," + String.valueOf(getBounds().getUpperLimit(0)) + "]");

        for (int i = 1; i < getDimension(); i++) {
            buffer.append(", " + String.valueOf(getParameterValue(i)));
            buffer.append("[" + String.valueOf(getBounds().getLowerLimit(i)));
            buffer.append("," + String.valueOf(getBounds().getUpperLimit(i)) + "]");
        }
        return buffer.toString();
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return COMPOUND_PARAMETER;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            CompoundParameter compoundParameter = new CompoundParameter(COMPOUND_PARAMETER);

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
            return CompoundParameter.class;
        }
    };

    // ****************************************************************
    // Parameter listener interface
    // ****************************************************************

    public void parameterChangedEvent(Parameter parameter, int index) {

        int dim = 0;
        for (int i = 0; i < parameters.length; i++) {
            if (parameter == parameters[i]) {
                fireParameterChangedEvent(dim + index);
                break;
            }
            dim += parameters[i].getDimension();
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


    private Parameter[] parameters = null;
    private IntersectionBounds bounds = null;
    private int dimension;
//	private String name;
}

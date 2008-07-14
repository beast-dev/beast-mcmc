/*
 * CompoundParameter.java
 *
 * Copyright (C) 2002-2007 Alexei Drummond and Andrew Rambaut
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

import java.util.ArrayList;
import java.util.List;

/**
 * A multidimensional parameter constructed from its component parameters.
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @version $Id: CompoundParameter.java,v 1.13 2005/06/14 10:40:34 rambaut Exp $
 */
public class CompoundParameter extends Parameter.Abstract implements ParameterListener {

    public static final String COMPOUND_PARAMETER = "compoundParameter";

    public CompoundParameter(Parameter[] params) {

        dimension = 0;
        for (Parameter parameter : parameters) {
            dimension += parameter.getDimension();
            parameter.addParameterListener(this);
        }

        this.parameters = new Parameter[dimension];
        this.pindex = new int[dimension];
        int k = 0;
        for (Parameter parameter : params) {
            for (int j = 0; j < parameter.getDimension(); j++) {
                parameters[j + k] = parameter;
                pindex[j + k] = j;
            }
            k += parameter.getDimension();
            uniqueParameters.add(parameter);
        }

    }

    public CompoundParameter(String name) {
        if (name != null) setId(name);
        dimension = 0;
    }

    public void addParameter(Parameter param) {

        uniqueParameters.add(param);
        if (parameters == null) {
            parameters = new Parameter[param.getDimension()];
            this.pindex = new int[param.getDimension()];
            for (int j = 0; j < param.getDimension(); j++) {
                parameters[j] = param;
                pindex[j] = j;
            }
        } else {
            Parameter[] newParams = new Parameter[parameters.length + param.getDimension()];
            int[] newIndices = new int[pindex.length + param.getDimension()];
            System.arraycopy(parameters, 0, newParams, 0, parameters.length);
            System.arraycopy(pindex, 0, newIndices, 0, pindex.length);

            for (int j = 0; j < param.getDimension(); j++) {
                newParams[j + parameters.length] = param;
                newIndices[j + pindex.length] = j;
            }

            parameters = newParams;
            pindex = newIndices;
        }
        dimension += param.getDimension();
        if (dimension != parameters.length) throw new RuntimeException();
        param.addParameterListener(this);
    }

    public final String getParameterName() {
        return getId();
    }

    public Parameter getParameter(int index) {
        return uniqueParameters.get(index);
    }

    public int getNumberOfParameters() {
        return uniqueParameters.size();
    }

    public final String getDimensionName(int dim) {

        return parameters[dim].getDimensionName(pindex[dim]);
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
        return parameters[dim].getParameterValue(pindex[dim]);
    }

    public double[] inspectParametersValues() {
        return getParameterValues();
    }

    public void setParameterValue(int dim, double value) {
        parameters[dim].setParameterValue(pindex[dim], value);
    }

    public void setParameterValueQuietly(int dim, double value) {
        parameters[dim].setParameterValueQuietly(pindex[dim], value);
    }

    protected void storeValues() {
        for (Parameter parameter : uniqueParameters) {
            parameter.storeParameterValues();
        }
    }

    protected void restoreValues() {
        for (Parameter parameter : uniqueParameters) {
            parameter.restoreParameterValues();
        }
    }

    protected final void acceptValues() {
        for (Parameter parameter : uniqueParameters) {
            parameter.acceptParameterValues();
        }
    }

    protected final void adoptValues(Parameter source) {
        // the parameters that make up a compound parameter will have
        // this function called on them individually so we don't need
        // to do anything here.
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

            CompoundParameter compoundParameter = new CompoundParameter((String) null);

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
        for (Parameter parameter1 : uniqueParameters) {
            if (parameter == parameter1) {
                fireParameterChangedEvent(dim + index);
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
            return parameters[dim].getBounds().getUpperLimit(pindex[dim]);
        }

        public double getLowerLimit(int dim) {
            return parameters[dim].getBounds().getLowerLimit(pindex[dim]);
        }

        public int getBoundsDimension() {
            return getDimension();
        }
    }


    private List<Parameter> uniqueParameters = new ArrayList<Parameter>();

    private Parameter[] parameters = null;
    private int[] pindex = null;
    private IntersectionBounds bounds = null;
    private int dimension;
}

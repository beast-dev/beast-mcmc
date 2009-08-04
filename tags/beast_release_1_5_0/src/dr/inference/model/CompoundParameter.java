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
public class CompoundParameter extends Parameter.Abstract implements VariableListener {

    public static final String COMPOUND_PARAMETER = "compoundParameter";

    public CompoundParameter(Parameter[] params) {

        dimension = 0;
        for (Parameter parameter : params) {
            dimension += parameter.getDimension();
            parameter.addParameterListener(this);
        }

        for (Parameter parameter : params) {
            for (int j = 0; j < parameter.getDimension(); j++) {
                parameters.add(parameter);
                pindex.add(j);
            }
            uniqueParameters.add(parameter);
        }

    }

    public CompoundParameter(String name) {
        this.name = name;
        dimension = 0;
    }

    public void addParameter(Parameter param) {

        uniqueParameters.add(param);
        for (int j = 0; j < param.getDimension(); j++) {
            parameters.add(param);
            pindex.add(j);
        }
        dimension += param.getDimension();
        if (dimension != parameters.size()) {
            throw new RuntimeException(
                    "dimension=" + dimension + " parameters.size()=" + parameters.size()
            );
        }
        param.addParameterListener(this);
    }

    public void removeParameter(Parameter param) {

        int dim = 0;
        for (Parameter parameter : uniqueParameters) {
            if (parameter == param) {
                break;
            }
            dim += parameter.getDimension();
        }

        for (int i = 0; i < param.getDimension(); i++) {
            parameters.remove(dim);
            pindex.remove(dim);
        }

        if (parameters.contains(param)) throw new RuntimeException();

        uniqueParameters.remove(param);

        dimension -= param.getDimension();
        if (dimension != parameters.size()) throw new RuntimeException();
        param.removeParameterListener(this);
    }

    /**
     * @return name if the parameter has been given a specific name, else it returns getId()
     */
    public final String getParameterName() {
        if (name != null) return name;
        return getId();
    }

    public Parameter getParameter(int index) {
        return uniqueParameters.get(index);
    }

    public int getNumberOfParameters() {
        return uniqueParameters.size();
    }

    public final String getDimensionName(int dim) {

        return parameters.get(dim).getDimensionName(pindex.get(dim));
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
        Parameter p = parameters.get(index);
        int pi = pindex.get(index);

        parameters.add(index, p);
        pindex.add(index, pi);

        p.addDimension(pi, value);
        for (int i = pi; i < p.getDimension(); i++) {
            pindex.set(index, i);
            index += 1;
        }
    }

    public double removeDimension(int index) {
        Parameter p = parameters.get(index);
        int pi = pindex.get(index);

        parameters.remove(index);
        pindex.remove(index);

        double v = p.removeDimension(pi);
        for (int i = pi; i < p.getDimension(); i++) {
            pindex.set(index, i);
            index += 1;
        }
        return v;
    }

    private void createBounds() {
        bounds = new IntersectionBounds(getDimension());
        bounds.addBounds(new CompoundBounds());
    }

    public double getParameterValue(int dim) {
        return parameters.get(dim).getParameterValue(pindex.get(dim));
    }

    public double[] inspectParametersValues() {
        return getParameterValues();
    }

    public void setParameterValue(int dim, double value) {
        parameters.get(dim).setParameterValue(pindex.get(dim), value);
    }

    public void setParameterValueQuietly(int dim, double value) {
        parameters.get(dim).setParameterValueQuietly(pindex.get(dim), value);
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

        private final XMLSyntaxRule[] rules;{
            rules = new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class, 1, Integer.MAX_VALUE),
            };
        }

        public Class getReturnType() {
            return CompoundParameter.class;
        }
    };

    // ****************************************************************
    // Parameter listener interface
    // ****************************************************************

    public void variableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {

        int dim = 0;
        for (Parameter parameter1 : uniqueParameters) {
            if (variable == parameter1) {
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

        public Double getUpperLimit(int dim) {
            return parameters.get(dim).getBounds().getUpperLimit(pindex.get(dim));
        }

        public Double getLowerLimit(int dim) {
            return parameters.get(dim).getBounds().getLowerLimit(pindex.get(dim));
        }

        public int getBoundsDimension() {
            return getDimension();
        }
    }

    private final List<Parameter> uniqueParameters = new ArrayList<Parameter>();

    private final ArrayList<Parameter> parameters = new ArrayList<Parameter>();
    private final ArrayList<Integer> pindex = new ArrayList<Integer>();
    private IntersectionBounds bounds = null;
    private int dimension;
    private String name;

    public static void main(String[] args) {

        Parameter param1 = new Parameter.Default(2);
        Parameter param2 = new Parameter.Default(2);
        Parameter param3 = new Parameter.Default(2);

        System.out.println(param1.getDimension());

        CompoundParameter parameter = new CompoundParameter(new Parameter[]{param1, param2});
        parameter.addParameter(param3);
        parameter.removeParameter(param2);
    }
}

/*
 * CompoundParameter.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

    public CompoundParameter(String name, Parameter[] params) {
        this(name);
        for (Parameter parameter : params) {
            dimension += parameter.getDimension();
            parameter.addParameterListener(this);
        }

        for (Parameter parameter : params) {
            for (int j = 0; j < parameter.getDimension(); j++) {
                parameters.add(parameter);
                pIndex.add(j);
            }
            uniqueParameters.add(parameter);
            labelParameter(parameter);
        }
    }

    public CompoundParameter(String name) {
        this.name = name;
        dimension = 0;
    }

    private void labelParameter(Parameter parameter) {
        if (parameter.getParameterName() == null) {
            String parameterName = name + uniqueParameters.size();
            parameter.setId(parameterName);
        }
    }

    public void addParameter(Parameter param) {

        uniqueParameters.add(param);
        for (int j = 0; j < param.getDimension(); j++) {
            parameters.add(param);
            pIndex.add(j);
        }
        dimension += param.getDimension();
        if (dimension != parameters.size()) {
            throw new RuntimeException(
                    "dimension=" + dimension + " parameters.size()=" + parameters.size()
            );
        }
        param.addParameterListener(this);
        labelParameter(param);
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
            pIndex.remove(dim);
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

    public int getParameterCount() {
        return uniqueParameters.size();
    }

    public String getDimensionName(int dim) {
        return parameters.get(dim).getDimensionName(pIndex.get(dim));
    }

    public int getDimension() {
        return dimension;
    }

    public void setDimension(int dim) {
        throw new UnsupportedOperationException();
    }

    public void addBounds(Bounds<Double> boundary) {

        if (bounds == null) {
            bounds = new CompoundBounds();
//            return;
        } //else {
        IntersectionBounds newBounds = new IntersectionBounds(getDimension());
        newBounds.addBounds(bounds);
        newBounds.addBounds(boundary);

//        }
//        ((IntersectionBounds) bounds).addBounds(boundary);
        bounds = newBounds;
    }

    public Bounds<Double> getBounds() {

        if (bounds == null) {
            bounds = new CompoundBounds();
        }
        return bounds;
    }

    public void addDimension(int index, double value) {
        Parameter p = parameters.get(index);
        int pi = pIndex.get(index);

        parameters.add(index, p);
        pIndex.add(index, pi);

        p.addDimension(pi, value);
        for (int i = pi; i < p.getDimension(); i++) {
            pIndex.set(index, i);
            index += 1;
        }
    }

    public double removeDimension(int index) {
        Parameter p = parameters.get(index);
        int pi = pIndex.get(index);

        parameters.remove(index);
        pIndex.remove(index);

        double v = p.removeDimension(pi);
        for (int i = pi; i < p.getDimension(); i++) {
            pIndex.set(index, i);
            index += 1;
        }
        return v;
    }

    public void fireParameterChangedEvent() {
        doNotPropagateChangeUp = true;
        for (Parameter p : parameters) {
            p.fireParameterChangedEvent();
        }
        doNotPropagateChangeUp = false;
        fireParameterChangedEvent(-1, ChangeType.ALL_VALUES_CHANGED);
    }

    public double getParameterValue(int dim) {
        return parameters.get(dim).getParameterValue(pIndex.get(dim));
    }

    public void setParameterValue(int dim, double value) {
        parameters.get(dim).setParameterValue(pIndex.get(dim), value);
    }

//    public void setParameterValue(int row, int column, double a)
//    {
//        getParameter(column).setParameterValue(row, a);
//    }

    public void setParameterValueQuietly(int dim, double value) {
        parameters.get(dim).setParameterValueQuietly(pIndex.get(dim), value);
    }

//    public void setParameterValueQuietly(int row, int column, double a){
//        getParameter(column).setParameterValueQuietly(row, a);
//    }

    public void setParameterValueNotifyChangedAll(int dim, double value) {
        parameters.get(dim).setParameterValueNotifyChangedAll(pIndex.get(dim), value);
    }

//    public void setParameterValueNotifyChangedAll(int row, int column, double val){
//        getParameter(column).setParameterValueNotifyChangedAll(row, val);
//    }

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

    protected void acceptValues() {
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
        return toStringCompoundParameter(getDimension());
    }

    protected String toStringCompoundParameter(int dim) {
        StringBuilder buffer = new StringBuilder(String.valueOf(getParameterValue(0)));
        final Bounds bounds = getBounds();
        buffer.append("[").append(String.valueOf(bounds.getLowerLimit(0)));
        buffer.append(",").append(String.valueOf(bounds.getUpperLimit(0))).append("]");

        for (int i = 1; i < dim; i++) {
            buffer.append(", ").append(String.valueOf(getParameterValue(i)));
            buffer.append("[").append(String.valueOf(bounds.getLowerLimit(i)));
            buffer.append(",").append(String.valueOf(bounds.getUpperLimit(i))).append("]");
        }
        return buffer.toString();
    }

    // ****************************************************************
    // Parameter listener interface
    // ****************************************************************

    public void variableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {

        int dim = 0;
        for (Parameter parameter1 : uniqueParameters) {
            if (variable == parameter1) {
                if (!doNotPropagateChangeUp) {
                    int subparameterIndex = (index == -1) ? -1 : dim + index;
                    fireParameterChangedEvent(subparameterIndex, type);
                }
                break;
            }
            dim += parameter1.getDimension();
        }
    }

    public double getParameterValue(int index, int parameter) {
        return getParameter(parameter).getParameterValue(index);
    }

    // ****************************************************************
    // Private and protected stuff
    // ****************************************************************

    private class CompoundBounds implements Bounds<Double> {

        public Double getUpperLimit(int dim) {
            return parameters.get(dim).getBounds().getUpperLimit(pIndex.get(dim));
        }

        public Double getLowerLimit(int dim) {
            return parameters.get(dim).getBounds().getLowerLimit(pIndex.get(dim));
        }

        public int getBoundsDimension() {
            return getDimension();
        }
    }

    protected ArrayList<Parameter> getParameters(){
        return parameters;
    }

    private final List<Parameter> uniqueParameters = new ArrayList<Parameter>();

    private final ArrayList<Parameter> parameters = new ArrayList<Parameter>();
    private final ArrayList<Integer> pIndex = new ArrayList<Integer>();
    private Bounds<Double> bounds = null;
    private int dimension;
    private String name;

    private boolean doNotPropagateChangeUp = false;

    public static void main(String[] args) {

        Parameter param1 = new Parameter.Default(2);
        Parameter param2 = new Parameter.Default(2);
        Parameter param3 = new Parameter.Default(2);

        System.out.println(param1.getDimension());

        CompoundParameter parameter = new CompoundParameter("parameter", new Parameter[]{param1, param2});
        parameter.addParameter(param3);
        parameter.removeParameter(param2);
    }
}

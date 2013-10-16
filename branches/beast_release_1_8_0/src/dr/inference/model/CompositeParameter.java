/*
 * CompositeParameter.java
 *
 * Copyright (C) 2002-2013 Alexei Drummond, Andrew Rambaut & Marc A. Suchard
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
 * A parameter which is the sum of a set of other parameters
 *
 * @author Andrew Rambaut
 * @version $Id$
 */
public class CompositeParameter extends Parameter.Abstract implements VariableListener {

    public CompositeParameter(String name, Parameter[] params) {
        this(name);

        for (Parameter parameter : params) {
            addParameter(parameter);
        }
    }

    public CompositeParameter(String name) {
        this.name = name;
        dimension = 0;
    }

    public void addParameter(Parameter param) {

        if (dimension == 0) {
            dimension = param.getDimension();
        }
        if (param.getDimension() != dimension && param.getDimension() != 1) {
            throw new RuntimeException("subsequent parameters do not match the dimensionality of the first (or 1)");
        }

        param.addParameterListener(this);
        parameters.add(param);
    }

    /**
     * @return name if the parameter has been given a specific name, else it returns getId()
     */
    public final String getParameterName() {
        if (name != null) return name;
        return getId();
    }

    public final String getDimensionName(int dim) {
        return parameters.get(0).getDimensionName(dim);
    }

    public int getDimension() {
        return dimension;
    }

    public void setDimension(int dim) {
        throw new RuntimeException();
    }

    public void addBounds(Bounds<Double> bounds) {
        throw new RuntimeException("Can't add bounds to a composite parameter, only its components");
    }

    public Bounds<Double> getBounds() {
        return null;
    }

    public void addDimension(int index, double value) {
        throw new RuntimeException();
    }

    public double removeDimension(int index) {
        throw new RuntimeException();
    }


    public double getParameterValue(int dim) {
        double value = parameters.get(0).getParameterValue(dim);
        for (int i = 1; i < parameters.size(); i++) {
            if (parameters.get(i).getDimension() == 1) {
                value += parameters.get(i).getParameterValue(0);
            } else {
                value += parameters.get(i).getParameterValue(dim);
            }
        }
        return value;
    }

    public void setParameterValue(int dim, double value) {
        throw new RuntimeException("Can't operate on a composite parameter, only its components");
    }

    public void setParameterValueQuietly(int dim, double value) {
        throw new RuntimeException("Can't operate on a composite parameter, only its components");
    }

    public void setParameterValueNotifyChangedAll(int dim, double value){
        throw new RuntimeException("Can't operate on a composite parameter, only its components");
    }

    protected void storeValues() {
        for (Parameter parameter : parameters) {
            parameter.storeParameterValues();
        }
    }

    protected void restoreValues() {
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

    // ****************************************************************
    // Parameter listener interface
    // ****************************************************************

    public void variableChangedEvent(Variable variable, int index, ChangeType type) {

        if (variable.getSize() > 1) {
            fireParameterChangedEvent(index, type);
        } else {
            fireParameterChangedEvent();
        }
    }

    private final ArrayList<Parameter> parameters = new ArrayList<Parameter>();

    private int dimension;
    private String name;

}

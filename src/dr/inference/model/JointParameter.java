/*
 * JointParameter.java
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

/**
 * A parameter which controls the values of a set of other parameters
 *
 * @author Andrew Rambaut
 * @version $Id$
 */
public class JointParameter extends Parameter.Abstract implements VariableListener {

    public JointParameter(String name, Parameter[] params) {
        this(name);

        for (Parameter parameter : params) {
            addParameter(parameter);
        }
    }

    public JointParameter(String name) {
        this.name = name;
        dimension = 0;
    }

    public void addParameter(Parameter param) {

        if (dimension == 0) {
            dimension = param.getDimension();
        } else {
            for (int dim = 0; dim < dimension; dim++) {
                param.setParameterValue(dim, parameters.get(0).getParameterValue(dim));
            }
            if (param.getDimension() != dimension) {
                throw new RuntimeException("subsequent parameters do not match the dimensionality of the first");
            }
        }

        // AR - I think we ignore the messages from the containing parameters. Possibly would be a good
        // idea to check they don't change independently of this one.
//        param.addParameterListener(this);

        parameters.add(param);
    }

    /**
     * @return name if the parameter has been given a specific name, else it returns getId()
     */
    public final String getParameterName() {
        if (name != null) return name;
        return getId();
    }

    public int getDimension() {
        return dimension;
    }

    public void setDimension(int dim) {
        throw new RuntimeException();
    }

    public void addBounds(Bounds<Double> bounds) {
        throw new RuntimeException("Can't add bounds to a joint parameter, only its components");
    }

    public Bounds<Double> getBounds() {
        return parameters.get(0).getBounds();
    }

    public void addDimension(int index, double value) {
        throw new RuntimeException();
    }

    public double removeDimension(int index) {
        throw new RuntimeException();
    }


    public double getParameterValue(int dim) {
        return parameters.get(0).getParameterValue(dim);
    }

    public void setParameterValue(int dim, double value) {
        for (Parameter parameter : parameters) {
            parameter.setParameterValue(dim, value);
        }
    }

    public void setParameterValueQuietly(int dim, double value) {
        for (Parameter parameter : parameters) {
            parameter.setParameterValueQuietly(dim, value);
        }
    }

    public void setParameterValueNotifyChangedAll(int dim, double value){
        for (Parameter parameter : parameters) {
            parameter.setParameterValueNotifyChangedAll(dim, value);
        }
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

/*
 * TransformedParameter.java
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

import dr.util.Transform;

/**
 * @author Marc A. Suchard
 * @author Paul Bastide
 */

public class TransformedParameter extends Parameter.Abstract implements VariableListener {

    public TransformedParameter(Parameter parameter, Transform transform) {
        this(parameter, transform, false);
    }

    public TransformedParameter(Parameter parameter, Transform transform, boolean inverse) {
        this.parameter = parameter;
        this.transform = transform;
        this.inverse = inverse;
        this.parameter.addVariableListener(this);
        Bounds bounds = parameter.getBounds();
        if (bounds != null && !(transform instanceof Transform.MultivariateTransform)) {
            addBounds(bounds);
        }
    }

    public int getDimension() {
        return parameter.getDimension();
    }

    protected void storeValues() {
        parameter.storeParameterValues();
    }

    protected void restoreValues() {
        parameter.restoreParameterValues();
    }

    protected void acceptValues() {
        parameter.acceptParameterValues();
    }

    protected void adoptValues(Parameter source) {
        parameter.adoptParameterValues(source);
    }

    private double transform(double value) {
        if (inverse) {
            return transform.inverse(value);
        } else {
            return transform.transform(value);
        }
    }

    protected double[] transform(double values[]) {
        if (inverse) {
            return transform.inverse(values, 0, values.length);
        } else {
            return transform.transform(values, 0, values.length);
        }
    }

    private double inverse(double value) {
        return inverse ? transform.transform(value) : transform.inverse(value);
    }

    protected double[] inverse(double values[]) {
        if (inverse) {
            return transform.transform(values, 0, values.length);
        } else {
            return transform.inverse(values, 0, values.length);
        }
    }

    public double getParameterValue(int dim) {
        return transform(parameter.getParameterValue(dim));
    }

    public double getParameterUntransformedValue(int dim) {
        return parameter.getParameterValue(dim);
    }

    public double[] getParameterUntransformedValues() {
        return parameter.getParameterValues();
    }

    public void setParameterUntransformedValue(int dim, double a) {
        parameter.setParameterValue(dim, a);
    }

    public Parameter getUntransformedParameter() {
        return parameter;
    }

    public void setParameterValue(int dim, double value) {
        parameter.setParameterValue(dim, inverse(value));
    }

    public void setParameterValueQuietly(int dim, double value) {
        parameter.setParameterValueQuietly(dim, inverse(value));
    }

    public void setParameterValueNotifyChangedAll(int dim, double value) {
        parameter.setParameterValueNotifyChangedAll(dim, inverse(value));
    }

    public String getParameterName() {
        if (getId() == null)
            return "transformed." + parameter.getParameterName();
        return getId();
    }

    public void addBounds(Bounds<Double> bounds) {
        final int dim = bounds.getBoundsDimension();
        final double[] lower = new double[dim];
        final double[] upper = new double[dim];
        for (int i = 0; i < dim; ++i) {
            final double transformedLowerBound = transform(bounds.getLowerLimit(i));
            final double transformedUpperBound = transform(bounds.getUpperLimit(i));
            if (transformedLowerBound < transformedUpperBound) {
                lower[i] = transformedLowerBound;
                upper[i] = transformedUpperBound;
            } else {
                lower[i] = transformedUpperBound;
                upper[i] = transformedLowerBound;
            }
        }
        transformedBounds = new DefaultBounds(upper, lower);

//        System.err.println("Started with:");
//        for (int i = 0; i < dim; ++i) {
//            System.err.print("\t" + bounds.getLowerLimit(i));
//        }
//        System.err.println("");
//        for (int i = 0; i < dim; ++i) {
//            System.err.print("\t" + bounds.getUpperLimit(i));
//        }
//        System.err.println("\n");
//
//        System.err.println("Ended with:");
//        for (int i = 0; i < dim; ++i) {
//            System.err.print("\t" + transformedBounds.getLowerLimit(i));
//        }
//        System.err.println("");
//        for (int i = 0; i < dim; ++i) {
//            System.err.print("\t" + transformedBounds.getUpperLimit(i));
//        }
//        System.err.println("\n");
//
//        parameter.addBounds(transformedBounds);
//        throw new RuntimeException("Should not call addBounds() on transformed parameter");
    }

    public boolean check() {
        return parameter.check();
    }

    public Bounds<Double> getBounds() {
        return transformedBounds;
//        throw new RuntimeException("Should not call getBounds() on transformed parameter");
    }

    public void addDimension(int index, double value) {
        throw new RuntimeException("Not yet implemented.");
    }

    public double removeDimension(int index) {
        throw new RuntimeException("Not yet implemented.");
    }

    @Override
    public void fireParameterChangedEvent() {

        doNotPropagateChangeUp = true;
        parameter.fireParameterChangedEvent();
        doNotPropagateChangeUp = false;

        fireParameterChangedEvent(-1, ChangeType.ALL_VALUES_CHANGED);
    }

    @Override
    public void variableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        if (!doNotPropagateChangeUp) {
            fireParameterChangedEvent(index, type);
        }
    }

    public double diffLogJacobian(double[] oldValues, double[] newValues) {
        // Takes **untransformed** values
        if (inverse) {
            return -transform.getLogJacobian(transform(oldValues), 0, oldValues.length)
                    + transform.getLogJacobian(transform(newValues), 0, newValues.length);
        } else {
            return transform.getLogJacobian(oldValues, 0, oldValues.length)
                    - transform.getLogJacobian(newValues, 0, newValues.length);
        }
    }

    public Transform getTransform() {
        return transform;
    }

    protected final Parameter parameter;
    protected final Transform transform;
    protected final boolean inverse;
    protected Bounds<Double> transformedBounds;

    protected boolean doNotPropagateChangeUp = false;
}
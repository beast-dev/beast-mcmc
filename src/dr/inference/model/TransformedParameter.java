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

    private double inverse(double value) {
        return inverse ? transform.transform(value) : transform.inverse(value);
    }

    public double getParameterValue(int dim) {
        return transform(parameter.getParameterValue(dim));
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
            lower[i] = inverse(bounds.getLowerLimit(i));
            upper[i] = inverse(bounds.getUpperLimit(i));
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
        parameter.addBounds(transformedBounds);
//        throw new RuntimeException("Should not call addBounds() on transformed parameter");
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

    public void variableChangedEvent(Variable variable, int index, ChangeType type) {
        // Propogate change up model graph
        fireParameterChangedEvent(index, type);
    }

    private final Parameter parameter;
    private final Transform transform;
    private final boolean inverse;
    private Bounds<Double> transformedBounds = null;
}

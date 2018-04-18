/*
 * TransformedMultivariateParameter.java
 *
 * Copyright (c) 2002-2018 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

public class TransformedMultivariateParameter extends TransformedParameter {

//    private double[] transformedValues;
//    private double[] storedTransformedValues; //TODO store/restore mechanism for TransformedParameter ?

    public TransformedMultivariateParameter(Parameter parameter, Transform.MultivariableTransform transform) {
        this(parameter, transform, false);
    }

    public TransformedMultivariateParameter(Parameter parameter, Transform.MultivariableTransform transform, boolean inverse) {
        super(parameter, transform, inverse);
//        this.transformedValues = transform(parameter.getParameterValues());
    }

    public double getParameterValue(int dim) {
        return transform(parameter.getParameterValues())[dim];
    }

    public void setParameterValue(int dim, double value) {
        double[] transformedValues = transform(parameter.getParameterValues());
        transformedValues[dim] = value;
        double[] newValues = inverse(transformedValues);
        // Need to update all values
        parameter.setParameterValueNotifyChangedAll(0, newValues[0]); // Warn everyone is changed
        for (int i = 1; i < parameter.getDimension(); i++) {
            parameter.setParameterValueQuietly(i, newValues[i]); // Do the rest quietly
        }
    }

    public void setParameterValueQuietly(int dim, double value) {
        double[] transformedValues = transform(parameter.getParameterValues());
        transformedValues[dim] = value;
        double[] newValues = inverse(transformedValues);
        // Need to update all values
        for (int i = 0; i < parameter.getDimension(); i++) {
            parameter.setParameterValueQuietly(i, newValues[i]);
        }
    }

    public void setParameterValueNotifyChangedAll(int dim, double value) {
        setParameterValue(dim, value);
    }

    public void addBounds(Bounds<Double> bounds) {
//        parameter.addBounds(new DefaultBounds(null, null));
        throw new RuntimeException("Should not call addBounds() on transformed parameter");
    }

    public Bounds<Double> getBounds() {
//        throw new RuntimeException("Should not call addBounds() on transformed parameter");
        return new DefaultBounds(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, parameter.getDimension());
    }

//    @Override
//    protected final void storeValues() {
//        System.err.println("YOUPI!");
//        super.storeValues();
//        if (storedTransformedValues == null || storedTransformedValues.length != transformedValues.length) {
//            storedTransformedValues = new double[transformedValues.length];
//        }
//        System.arraycopy(transformedValues, 0, storedTransformedValues, 0, storedTransformedValues.length);
//    }
//
//    @Override
//    protected final void restoreValues() {
//        System.err.println("YOUPI!");
//        super.restoreValues();
//        //swap the arrays
//        double[] temp = storedTransformedValues;
//        storedTransformedValues = transformedValues;
//        transformedValues = temp;
//    }
}

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

    private double[] transformedValues;
    private double[] unTransformedValues;
//    private double[] storedTransformedValues; //TODO store/restore mechanism for TransformedParameter ?

    private boolean valuesKnown = false;

    public TransformedMultivariateParameter(Parameter parameter, Transform.MultivariableTransform transform) {
        this(parameter, transform, false);
    }

    public TransformedMultivariateParameter(Parameter parameter, Transform.MultivariableTransform transform, boolean inverse) {
        super(parameter, transform, inverse);
        this.unTransformedValues = parameter.getParameterValues();
        this.transformedValues = transform(unTransformedValues);
    }

    public double getParameterValue(int dim) {
        update();
        return transformedValues[dim];
    }

    protected void storeValues() {
        super.storeValues();
    }

    protected void restoreValues() {
        super.restoreValues();
        valuesKnown = false;
    }

//    public void variableChangedEvent(Variable variable, int index, ChangeType type) {
//        valuesKnown = false;
//        // Propogate change up model graph
//        fireParameterChangedEvent(index, type);
//    }

    public void setParameterValue(int dim, double value) {
        setParameterValueQuietly(dim, value);
        parameter.fireParameterChangedEvent();
    }

    public void setParameterValueQuietly(int dim, double value) {
        update();
        transformedValues[dim] = value;
        updateParameterQuietlyFromTransformedValues();
    }

    @Override
    public void setAllParameterValuesQuietly(double[] values) {
        if (values.length != transformedValues.length) {
            throw new IllegalArgumentException("supplied values must be of same dimension as transformed parameter");
        }

        for (int i = 0; i < transformedValues.length; i++) {
            transformedValues[i] = values[i];
        }
        updateParameterQuietlyFromTransformedValues();
    }

    private void updateParameterQuietlyFromTransformedValues() {
        unTransformedValues = inverse(transformedValues);
        // Need to update all values
        for (int i = 0; i < parameter.getDimension(); i++) {
            parameter.setParameterValueQuietly(i, unTransformedValues[i]);
        }
    }

    public void setParameterValueNotifyChangedAll(int dim, double value) {
        setParameterValue(dim, value);
    }

    public void addBounds(Bounds<Double> bounds) {
        // TODO: Check bounds of the parameter ?  XJ: bounds can be quite arbitrary in this case.  I decided to allow manual setup through parser.
        transformedBounds = bounds;
    }

    private void update() {
        if (!valuesKnown) {
            unTransformedValues = parameter.getParameterValues();
            transformedValues = transform(unTransformedValues);
            valuesKnown = true;
        }
    }

    private boolean hasChanged() {
        for (int i = 0; i < unTransformedValues.length; i++) {
            if (parameter.getParameterValue(i) != unTransformedValues[i]) {
                return true;
            }
        }
        return false;
    }

    @Override
    public void variableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        valuesKnown = false;
        if (!doNotPropagateChangeUp) {
            fireParameterChangedEvent(-1, ChangeType.ALL_VALUES_CHANGED); //if one dimension of the untransformed parameter changes, it is very likely that many dimensions of the transformed parameter change
        }
    }
}

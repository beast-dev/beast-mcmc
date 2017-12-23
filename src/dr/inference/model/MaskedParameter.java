/*
 * MaskedParameter.java
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

import java.util.Arrays;

/**
 * @author Marc A. Suchard
 */

public class MaskedParameter extends Parameter.Abstract implements VariableListener {

    public enum Signaling {
        NORMAL, NO_DEPENDENT;
    }

    private final Signaling signaling;

    public MaskedParameter(Parameter parameter, Parameter maskParameter, boolean ones) {
        this(parameter, maskParameter, ones, Signaling.NORMAL);
    }

    public MaskedParameter(Parameter parameter, Parameter maskParameter, boolean ones,
                           Signaling signaling) {
        this(parameter, signaling);
        addMask(maskParameter, ones);
    }

    public MaskedParameter(Parameter parameter) {
        this(parameter, Signaling.NORMAL);
    }

    public MaskedParameter(Parameter parameter, Signaling signaling) {
        this.parameter = parameter;
        parameter.addParameterListener(this);

        this.map = new int[parameter.getDimension()];
        this.inverseMap = new int[parameter.getDimension()];

        this.storedMap = new int[parameter.getDimension()];
        this.storedInverseMap = new int[parameter.getDimension()];

        for (int i = 0; i < map.length; i++) {
            map[i] = i;
            inverseMap[i] = i;
        }
        length = map.length;

        this.signaling = signaling;
    }

    public void addMask(Parameter maskParameter, boolean ones) {
        if (maskParameter.getDimension() != parameter.getDimension())
            throw new RuntimeException("Masking parameter '" + maskParameter.getId() + "' dimension must equal base parameter '" +
                    parameter.getId() + "' dimension");
        this.maskParameter = maskParameter;
        maskParameter.addParameterListener(this);
        if (ones)
            equalValue = 1;
        else
            equalValue = 0;
        updateMask();
    }

    private void updateMask() {
        length = updateMask(maskParameter, map, inverseMap, equalValue);
        bounds = null;
    }

    public static int updateMask(Parameter maskParameter, int[] map, int[] inverseMap, int equalValue) {
        int index = 0;
        for (int i = 0; i < maskParameter.getDimension(); i++) {
            // TODO Add a threshold attribute for continuous value masking
            final int maskValue = (int) maskParameter.getParameterValue(i);
            if (maskValue == equalValue) {
                map[index] = i;
                inverseMap[i] = index;
                index++;
            } else {
                inverseMap[i] = -1; // Keep track of indices from parameter than do NOT correspond to entries in mask
            }
        }
        return index;
    }

    public int getDimension() {
//        if (length == 0)
//            throw new RuntimeException("Zero-dimensional parameter!");
//        // TODO Need non-fatal mechanism to check for zero-dimensional parameters
        return length;
    }

    protected void storeValues() {
        parameter.storeParameterValues();
        maskParameter.storeParameterValues();

        System.arraycopy(map, 0, storedMap, 0, map.length);
        System.arraycopy(inverseMap, 0, storedInverseMap, 0, inverseMap.length);
    }

    protected void restoreValues() {
        parameter.restoreParameterValues();
        maskParameter.restoreParameterValues();

        int[] tmp = storedMap;
        storedMap = map;
        map = tmp;

        tmp = storedInverseMap;
        storedInverseMap = inverseMap;
        inverseMap = tmp;
    }

    public void fireParameterChangedEvent() {
        if (signaling == Signaling.NORMAL) {
            doNotPropagateChangeUp = true;
            parameter.fireParameterChangedEvent();
            doNotPropagateChangeUp = false;
        }
        super.fireParameterChangedEvent();
    }

    public void fireParameterChangedEvent(int index, Parameter.ChangeType type) {
        if (signaling == Signaling.NORMAL) {
            doNotPropagateChangeUp = true;
            parameter.fireParameterChangedEvent(index, type);
            doNotPropagateChangeUp = false;
        }
        super.fireParameterChangedEvent(index, type);
    }

    protected void acceptValues() {
        parameter.acceptParameterValues();
        maskParameter.acceptParameterValues();
    }

    protected void adoptValues(Parameter source) { throw new IllegalArgumentException("Not yet implemented"); }

    public double getParameterValue(int dim) {
//        if (!isMapValid) updateMask();
        return parameter.getParameterValue(map[dim]);
    }

    public void setParameterValue(int dim, double value) {
//        if (!isMapValid) updateMask();
        parameter.setParameterValue(map[dim], value);
    }

    public void setParameterValueQuietly(int dim, double value) {
        parameter.setParameterValueQuietly(map[dim], value);
    }

    public void setParameterValueNotifyChangedAll(int dim, double value) {
        parameter.setParameterValueNotifyChangedAll(map[dim], value);
    }

    @SuppressWarnings("unused")
    public double getParameterMaskValue(int i){
        return maskParameter.getParameterValue(i);
    }

    public String getParameterName() {
        if (getId() == null)
            return "masked." + parameter.getParameterName();
        return getId();
    }

    public void addBounds(Bounds<Double> bounds) {
        final int dimNotMasked = parameter.getDimension();
        final double[] lower = new double[dimNotMasked];
        final double[] upper = new double[dimNotMasked];

        Arrays.fill(lower, Double.NEGATIVE_INFINITY);
        Arrays.fill(upper, Double.POSITIVE_INFINITY);

        assert (bounds.getBoundsDimension() == getDimension());

        final int dimMasked = getDimension();
        for (int i = 0; i < dimMasked; ++i) {
            lower[map[i]] = bounds.getLowerLimit(i);
            upper[map[i]] = bounds.getUpperLimit(i);
        }

        DefaultBounds notMaskedBounds = new DefaultBounds(upper, lower);
        parameter.addBounds(notMaskedBounds);
    }

    public Bounds<Double> getBounds() {
        if (bounds == null) {

            // Create masked bounds
            Bounds<Double> oldBounds = parameter.getBounds();
            if (oldBounds != null) {
                double[] upper = new double[length];
                double[] lower = new double[length];

                for (int i = 0; i < length; ++i) {
                    upper[i] = oldBounds.getUpperLimit(map[i]);
                    lower[i] = oldBounds.getLowerLimit(map[i]);
                }

                bounds = new DefaultBounds(upper, lower);
            }
        }
        return bounds;
    }

    public Parameter getUnmaskedParameter() {
        return parameter;
    }

    public void addDimension(int index, double value) {
        throw new RuntimeException("Not yet implemented.");
    }

    public double removeDimension(int index) {
        throw new RuntimeException("Not yet implemented.");
    }

    public void variableChangedEvent(Variable variable, int index, ChangeType type) {
        if (variable == maskParameter) {
            updateMask();
            super.fireParameterChangedEvent();
        } else if (variable == parameter) { // variable == parameter
            if (!doNotPropagateChangeUp) {
                if (index == -1) {
                    super.fireParameterChangedEvent();
                } else if (inverseMap[index] != -1) {
                    super.fireParameterChangedEvent(inverseMap[index], type);
                }
            }
        } else {
            throw new IllegalArgumentException("Unknown variable");
        }
    }

    private final Parameter parameter;
    private Parameter maskParameter;

    private Bounds<Double> bounds = null;

    private int[] map;
    private int[] inverseMap;

    private int[] storedMap;
    private int[] storedInverseMap;

    private int length;
    private int equalValue;

    private boolean doNotPropagateChangeUp = false;
}

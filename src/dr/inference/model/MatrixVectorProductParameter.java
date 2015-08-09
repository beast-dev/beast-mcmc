/*
 * MatrixVectorProductParameter.java
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

/**
 * @author Marc Suchard
 */
public class MatrixVectorProductParameter extends Parameter.Abstract implements VariableListener {

    public MatrixVectorProductParameter(MatrixParameter matrix, Parameter vector) {
        // TODO Check dimensions
        this.matrix = matrix;
        this.vector = vector;

        matrix.addVariableListener(this);
        vector.addVariableListener(this);
    }

    public int getDimension() {
        return matrix.getRowDimension();
    }

    public void setDimension(int dim) {
        throwError("setDimension()");
    }

    protected void storeValues() {
        matrix.storeParameterValues();
        vector.storeParameterValues();
    }

    protected void restoreValues() {
        matrix.restoreParameterValues();
        vector.restoreVariableValues();
    }

    protected void acceptValues() {
        matrix.acceptParameterValues();
        vector.acceptParameterValues();
    }

    protected void adoptValues(Parameter source) {
        throw new RuntimeException("Not implemented");
    }

    public double getParameterValue(int dim) {
        double value = 0.0;
        for (int k = 0; k < matrix.getColumnDimension(); ++k) {
            value += matrix.getParameterValue(dim, k) * vector.getParameterValue(k);
        }
        return value;
    }

    public void setParameterValue(int dim, double value) {
        throwError("setParameterValue()");
    }

    public void setParameterValueQuietly(int dim, double value) {
        throwError("setParameterValueQuietly()");
    }

    public void setParameterValueNotifyChangedAll(int dim, double value) {
        throwError("setParameterValueNotifyChangedAll()");
    }

    private void throwError(String functionName) throws RuntimeException {
        throw new RuntimeException("Object " + getId() + " is a deterministic function. Calling "
                + functionName + " is not allowed");
    }

    public String getParameterName() {
        if (getId() == null) {
            StringBuilder sb = new StringBuilder("product");
            sb.append(".").append(matrix.getId());
            sb.append(".").append(vector.getId());
            setId(sb.toString());
        }
        return getId();
    }

    public void addBounds(Bounds bounds) {
        this.bounds = bounds;
    }

    public Bounds<Double> getBounds() {
        throwError("getBounds()");
        return null;
    }

    public void addDimension(int index, double value) {
        throw new RuntimeException("Not yet implemented.");
    }

    public double removeDimension(int index) {
        throw new RuntimeException("Not yet implemented.");
    }

    public void variableChangedEvent(Variable variable, int index, ChangeType type) {
        fireParameterChangedEvent(); // All dimensions may have changed
    }

    private final MatrixParameter matrix;
    private final Parameter vector;

    private Bounds bounds = null;
}

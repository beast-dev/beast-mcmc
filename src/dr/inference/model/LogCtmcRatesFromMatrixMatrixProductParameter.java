/*
 * LogCtmcRatesFromMatrixMatrixProductParameter.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.inference.model;

/**
 * @author Marc Suchard
 */
public class LogCtmcRatesFromMatrixMatrixProductParameter extends Parameter.Abstract implements VariableListener {

    public LogCtmcRatesFromMatrixMatrixProductParameter(MatrixParameter lhs,  MatrixParameter rhs) {

        assert lhs.getRowDimension() == rhs.getColumnDimension();
        assert rhs.getColumnDimension() == rhs.getRowDimension();

        this.lhs = lhs;
        this.rhs = rhs;
        this.dim = lhs.getRowDimension();

        this.map = makeMap(dim);

        lhs.addVariableListener(this);
        rhs.addVariableListener(this);

        Parameter.CONNECTED_PARAMETER_SET.add(lhs);
        Parameter.CONNECTED_PARAMETER_SET.add(rhs);
    }

    public int getDimension() {
        return dim;
    }

    public void setDimension(int dim) {
        throwError("setDimension()");
    }

    protected void storeValues() {
        lhs.storeParameterValues();
        rhs.storeParameterValues();
    }

    protected void restoreValues() {
        lhs.restoreParameterValues();
        rhs.restoreVariableValues();
    }

    protected void acceptValues() {
        lhs.acceptParameterValues();
        rhs.acceptParameterValues();
    }

    protected void adoptValues(Parameter source) {
        throwError("adoptValues");
    }

    public double getParameterValue(int dim) {

        final int row = map[dim][0];
        final int col = map[dim][1];

        double value = 0.0;
        for (int k = 0; k < lhs.getColumnDimension(); ++k) {
            value += lhs.getParameterValue(row, k) * rhs.getParameterValue(k, col);
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

    public boolean isImmutable() {
        return true;
    }

    private void throwError(String functionName) throws RuntimeException {
        throw new RuntimeException("Object " + getId() + " is a deterministic function. Calling "
                + functionName + " is not allowed");
    }

    public String getParameterName() {
        if (getId() == null) {
            String sb = "product" + "." + lhs.getId() +
                    "." + rhs.getId();
            setId(sb);
        }
        return getId();
    }

    public void addBounds(Bounds<Double> bounds) {
        this.bounds = bounds;
    }

    public Bounds<Double> getBounds() {
        return bounds;
    }

    public void addDimension(int index, double value) {
        throwError("addDimension");
    }

    public double removeDimension(int index) {
        throwError("removeDimension");
        return Double.NaN;
    }

    public void variableChangedEvent(Variable variable, int index, ChangeType type) {
        fireParameterChangedEvent(); // All dimensions may have changed
    }

    private int[][] makeMap(int ctmcDim) {
        int[][] map = new int[ctmcDim * (ctmcDim - 1) / 2][];

        int index = 0;
        for (int i = 1; i < ctmcDim; ++i) {
            for (int j = i + 2; j < ctmcDim; ++j) {
                map[index++] = new int[] {i, j};
            }
        }

        for (int j = 1; j < ctmcDim; ++j) {
            for (int i = j + 1; i < ctmcDim; ++i) {
                map[index++] = new int[] {i, j};
            }
        }

        return map;
    }

    private final MatrixParameter lhs;
    private final MatrixParameter rhs;
    private final int dim;
    private final int[][] map;

    private Bounds<Double> bounds;
}

/*
 * MatrixVectorProductParameter.java
 *
 * Copyright (c) 2002-2013 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
//        paramList = new ArrayList<Parameter>();
//        paramList.add(matrix);
//        paramList.add(vector);
        this.matrix = matrix;
        this.vector = vector;

//        for (Parameter p : paramList) {
//            p.addVariableListener(this);
//        }
        matrix.addVariableListener(this);
        vector.addVariableListener(this);
    }

    public int getDimension() {
//        return paramList.get(0).getDimension();
        return matrix.getRowDimension();
    }

    protected void storeValues() {
//        for (Parameter p : paramList) {
//            p.storeParameterValues();
//        }
        matrix.storeParameterValues();
        vector.storeParameterValues();
    }

    protected void restoreValues() {
//        for (Parameter p : paramList) {
//            p.restoreParameterValues();
//        }
        matrix.restoreParameterValues();
        vector.restoreVariableValues();
    }

    protected void acceptValues() {
//        for (Parameter p : paramList) {
//            p.acceptParameterValues();
//        }
        matrix.acceptParameterValues();
        vector.acceptParameterValues();
    }

    protected void adoptValues(Parameter source) {
        throw new RuntimeException("Not implemented");
    }

    public double getParameterValue(int dim) {
        double value = 0.0;
//        for (Parameter p : paramList) {
//            value *= p.getParameterValue(dim);
//        }
        for (int k = 0; k < matrix.getColumnDimension(); ++k) {
            value += matrix.getParameterValue(dim, k) * vector.getParameterValue(k);
        }
        return value;
    }

    public void setParameterValue(int dim, double value) {
        throw new RuntimeException("Not implemented");
    }

    public void setParameterValueQuietly(int dim, double value) {
        throw new RuntimeException("Not implemented");
    }

    public void setParameterValueNotifyChangedAll(int dim, double value) {
        throw new RuntimeException("Not implemented");
    }

    public String getParameterName() {
        if (getId() == null) {
            StringBuilder sb = new StringBuilder("product");
//            for (Parameter p : paramList) {
//                sb.append(".").append(p.getId());
//            }
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
        if (bounds == null) {
//            return paramList.get(0).getBounds(); // TODO
            throw new NullPointerException(getParameterName() + " parameter: Bounds not set");
//            return vector.getBounds();
        } else {
            return bounds;
        }
    }

    public void addDimension(int index, double value) {
        throw new RuntimeException("Not yet implemented.");
    }

    public double removeDimension(int index) {
        throw new RuntimeException("Not yet implemented.");
    }

    public void variableChangedEvent(Variable variable, int index, ChangeType type) {
        fireParameterChangedEvent(index, type);
    }

    //    private final List<Parameter> paramList;
    private final MatrixParameter matrix;
    private final Parameter vector;

    private Bounds bounds = null;
}

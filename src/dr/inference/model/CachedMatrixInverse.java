/*
 * CompoundSymmetricMatrix.java
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

import dr.math.matrixAlgebra.Matrix;
import dr.math.matrixAlgebra.ReadableMatrix;
import dr.math.matrixAlgebra.WrappedMatrix;
import dr.math.matrixAlgebra.WritableMatrix;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

/**
 * @author Marc Suchard
 * @author Paul Bastide
 */
public class CachedMatrixInverse extends CompoundParameter implements MatrixParameterInterface {

    private final MatrixParameterInterface base;
    private final int dim;

    private boolean inverseKnown;
    private boolean savedInverseKnown;

    private WrappedMatrix inverse;
    private WrappedMatrix savedInverse;

    public CachedMatrixInverse(String name,
                               MatrixParameterInterface base) {
        super(name, new Parameter[] { base });

        assert base.getColumnDimension() == base.getRowDimension();

        this.base = base;
        this.dim = base.getColumnDimension();
    }

    @Override
    public int getDimension() { return dim * dim; }

    @Override
    public double getParameterValue(int row, int col) {
        checkInverse();
        return inverse.get(row, col);
    }

    @Override
    public double getParameterValue(int index) {
        return getParameterValue(index / dim, index % dim);
    }

    @Override
    public double[][] getParameterAsMatrix() {

        double[][] matrix = new double[dim][dim];
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                matrix[i][j] = getParameterValue(i, j);
            }
        }

        return matrix;
    }

    @Override
    public void variableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        if (variable == base) {
            inverseKnown = false;
            fireParameterChangedEvent(index, type);
        } else {
            throw new IllegalArgumentException("Unknown variable");
        }
    }

    private static final boolean EMJL = true;

    private void computeInverse() {

        if (DEBUG) {
            System.err.println("CachedMatrixInverse.computeInverse()");
        }
        
        if (EMJL) {
            // TODO Avoid multiple copies
            DenseMatrix64F source = new DenseMatrix64F(base.getParameterAsMatrix());
            DenseMatrix64F destination = new DenseMatrix64F(getColumnDimension(), getColumnDimension());
            CommonOps.invert(source, destination);
            inverse = new WrappedMatrix.WrappedDenseMatrix(destination);
        } else {
         inverse = new WrappedMatrix.ArrayOfArray(new Matrix(base.getParameterAsMatrix()).inverse().toComponents());
        }
    }

    private void checkInverse() {
        if (!inverseKnown) {
            computeInverse();
            inverseKnown = true;
        }
    }

    @Override
    public int getColumnDimension() { return dim; }

    @Override
    public int getRowDimension() { return dim; }

    @Override
    public int getUniqueParameterCount() { return 1; }

    @Override
    public Parameter getUniqueParameter(int index) { return base; }

    @Override
    public void copyParameterValues(double[] destination, int offset) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public String toSymmetricString() {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public void setParameterValue(int index, double a) {
        throw new RuntimeException("Do not set entries directly");
    }

    @Override
    public void setParameterValue(int row, int column, double a) {
        throw new RuntimeException("Do not set entries directly");
    }

    @Override
    public void setParameterValueQuietly(int row, int col, double value) {
        throw new RuntimeException("Do not set entries directly");
    }

    @Override
    public void setParameterValueNotifyChangedAll(int row, int col, double value) {
        throw new RuntimeException("Do not set entries directly");
    }

    @Override
    public void setAllParameterValuesQuietly(double[] values, int offset) {
        throw new RuntimeException("Do not set entries directly");
    }

    @Override
    public double[] getColumnValues(int col) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public String getDimensionName(int index) {
        int dim = getColumnDimension();
        String row = Integer.toString(index / dim);
        String col = Integer.toString(index % dim);

        return getId() + row + col;
    }

    @Override
    protected void storeValues() {
        super.storeValues();
        savedInverseKnown = inverseKnown;

        if (inverse != null) {

            if (savedInverse == null) {
                savedInverse = new WrappedMatrix.WrappedDenseMatrix(new DenseMatrix64F(dim ,dim));
            }
            copy(inverse, savedInverse);
        }
    }

    @Override
    protected void restoreValues() {
        super.restoreValues();
        inverseKnown = savedInverseKnown;

        WrappedMatrix tmp = inverse;
        inverse = savedInverse;
        savedInverse = tmp;
    }

    static void copy(final ReadableMatrix source, final WritableMatrix destination) {
        assert source.getDim() == destination.getDim();

        final int length = source.getDim();
        for (int i = 0; i < length; ++i) {
            destination.set(i, source.get(i));
        }
    }

    @Override
    public String getReport() {
        return new WrappedMatrix.ArrayOfArray(getParameterAsMatrix()).toString();
    }

    private static final boolean DEBUG = true;

    public MatrixParameterInterface getBaseParameter() {
        return base;
    }
}

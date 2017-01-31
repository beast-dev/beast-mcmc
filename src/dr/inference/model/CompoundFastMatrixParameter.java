/*
 * CompoundFastMatrixParameter.java
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

import dr.xml.*;
import mpi.MPI;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Marc A. Suchard
 */
public class CompoundFastMatrixParameter extends CompoundParameter implements MatrixParameterInterface {

    private final int rowDimension;
    private int colDimension;

    private final List<MatrixParameterInterface> columns = new ArrayList<MatrixParameterInterface>();
    private final List<Integer> offsets = new ArrayList<Integer>();
    private final List<MatrixParameterInterface> matrices;

    public CompoundFastMatrixParameter(String name, List<MatrixParameterInterface> matrices) {
        super(name, compoundMatrices(matrices));

        this.matrices = matrices;

        rowDimension = matrices.get(0).getRowDimension();
        colDimension = 0;

        for (MatrixParameterInterface matrix : matrices) {
            if (matrix.getRowDimension() != rowDimension) {
                throw new IllegalArgumentException("Inconsistent row dimensions");
            }

            for (int i = 0; i < matrix.getColumnDimension(); ++i) {
                columns.add(matrix);
                offsets.add(i);
            }

            colDimension += matrix.getColumnDimension();
        }
    }

    private static Parameter[] compoundMatrices(List<MatrixParameterInterface> matrices) {
        int length = 0;
        for (MatrixParameterInterface matrix : matrices) {
            length += matrix.getUniqueParameterCount();
        }

        Parameter[] parameters = new Parameter[length];
        int index = 0;

        for (MatrixParameterInterface matrix : matrices) {
            for (int i = 0; i < matrix.getUniqueParameterCount(); ++i) {
                parameters[index] = matrix.getUniqueParameter(i);
                ++index;
            }
        }

        return parameters;
    }

    @Override
    public Parameter getParameter(int column) {
        return columns.get(column).getParameter(offsets.get(column));
    }

    @Override
    public double getParameterValue(int row, int col) {
        return columns.get(col).getParameterValue(row, offsets.get(col));
    }

    @Override
    public void setParameterValue(int row, int col, double value) {
        columns.get(col).setParameterValue(row, offsets.get(col), value);
    }

    @Override
    public void setParameterValueQuietly(int row, int col, double value) {
        columns.get(col).setParameterValueQuietly(row, offsets.get(col), value);
    }

    @Override
    public void setParameterValueNotifyChangedAll(int row, int column, double value) {
        columns.get(column).setParameterValueNotifyChangedAll(row, offsets.get(column), value);
    }

    @Override
    public double[] getColumnValues(int col) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public double[][] getParameterAsMatrix() {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public int getColumnDimension() {
        return colDimension;
    }

    @Override
    public int getRowDimension() {
        return rowDimension;
    }

    @Override
    public int getUniqueParameterCount() {
        return getParameterCount();
    }

    @Override
    public Parameter getUniqueParameter(int index) {
        return getParameter(index);
    }

    @Override
    public void copyParameterValues(double[] destination, int offset) {
        for (MatrixParameterInterface matrix : matrices) {
            matrix.copyParameterValues(destination, offset);
            offset += matrix.getRowDimension() * matrix.getColumnDimension();
        }
    }

    @Override
    public void setAllParameterValuesQuietly(double[] values, int offset) {
        for (MatrixParameterInterface matrix : matrices) {
            matrix.setAllParameterValuesQuietly(values, offset);
            offset += matrix.getRowDimension() * matrix.getColumnDimension();
        }
    }

    @Override
    public double[] getParameterValues() {
        int length = 0;
        for (MatrixParameterInterface matrix : matrices) {
            length += matrix.getRowDimension() * matrix.getColumnDimension();
        }
        double[] rtn = new double[length];
        copyParameterValues(rtn, 0);
        return rtn;
    }

    @Override
    public String toSymmetricString() {
        return MatrixParameter.toSymmetricString(this);
    }

    public final static String COMPOUND_FAST_MATRIX_PARAMETER = "compoundFastMatrixParameter";

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return COMPOUND_FAST_MATRIX_PARAMETER;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            List<MatrixParameterInterface> matrices = new ArrayList<MatrixParameterInterface>();

            for (int i = 0; i < xo.getChildCount(); ++i) {
                matrices.add((MatrixParameterInterface) xo.getChild(i));
            }

            final String name = xo.hasId() ? xo.getId() : null;

            return new CompoundFastMatrixParameter(name, matrices);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "A compound matrix parameter constructed from its component parameters.";
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                new ElementRule(MatrixParameterInterface.class, 1, Integer.MAX_VALUE),
        };

        public Class getReturnType() {
            return CompoundFastMatrixParameter.class;
        }
    };
}

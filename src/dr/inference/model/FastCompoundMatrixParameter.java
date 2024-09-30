/*
 * FastCompoundMatrixParameter.java
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

import dr.xml.*;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * @author Marc A. Suchard
 * @author Andrew Holbrokk
 */

public class FastCompoundMatrixParameter extends CompoundParameter implements MatrixParameterInterface {

    private final List<MatrixParameterInterface> matrices;
    private final int nRows;
    private final int nColumns;
    private final int[] matrixIndex;
    private final int[] columnIndex;
    private boolean doPropagateChangeUp = true;

    private FastCompoundMatrixParameter(String name, List<MatrixParameterInterface> matrices) {
        super(name);

        this.matrices = matrices;

        this.nRows = matrices.get(0).getRowDimension();
        int columns = 0;

        for (MatrixParameterInterface matrix : matrices) {
            if (matrix.getRowDimension() != nRows) {
                throw new RuntimeException("Invalid row dimensions");
            }

            columns += matrix.getColumnDimension();
        }

        this.nColumns = columns;
        this.matrixIndex = new int[nColumns];
        this.columnIndex = new int[nColumns];

        int offset = 0;
        for (int m = 0; m < matrices.size(); ++m) {
            MatrixParameterInterface matrix = matrices.get(m);
            int length = matrix.getColumnDimension();
            Arrays.fill(matrixIndex, offset, offset + length, m);
            for (int i = 0; i < length; ++i) {
                columnIndex[offset + i] = i;
            }
            offset += length;
        }


        for (MatrixParameterInterface matrix : matrices) {
            addParameter(matrix);
        }

    }

    @Override
    public void fireParameterChangedEvent() {
        fireParameterChangedEvent(-1, ChangeType.ALL_VALUES_CHANGED);
    }

    @Override
    public void fireParameterChangedEvent(int index, Parameter.ChangeType type) {

        assert (index == -1);
        assert (type == ChangeType.ALL_VALUES_CHANGED);

        doPropagateChangeUp = false;
        for (MatrixParameterInterface matrix : matrices) {
            matrix.fireParameterChangedEvent();
        }
        doPropagateChangeUp = true;
        super.fireParameterChangedEvent(index, type);
    }

    @Override
    public void variableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {

        assert (index == -1);
        assert (type == ChangeType.ALL_VALUES_CHANGED);

        if (doPropagateChangeUp) {
            fireParameterChangedEvent(-1, ChangeType.ALL_VALUES_CHANGED);
        }
    }

    @Override
    public Parameter getParameter(int index) {
        return matrices.get(matrixIndex[index]).getParameter(columnIndex[index]);
    }

//    @Override
//    public void setParameterValue(int index, double value) {
//        throw new RuntimeException("Not yet implemented");
//    }
//
//    @Override void setParameter

    @Override
    public void setParameterValue(int row, int col, double value) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public void setParameterValueQuietly(int row, int col, double value) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public void setParameterValueNotifyChangedAll(int row, int col, double value) {
        throw new RuntimeException("Not yet implemented");
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
        return nColumns;
    }

    @Override
    public int getRowDimension() {
        return nRows;
    }

    @Override
    public int getUniqueParameterCount() {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public Parameter getUniqueParameter(int index) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public void copyParameterValues(double[] destination, int offset) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public void setAllParameterValuesQuietly(double[] values, int offset) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public String toSymmetricString() {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public boolean isConstrainedSymmetric() {
        throw new RuntimeException("Not yet implemented");
    }

    private final static String FAST_COMPOUND_MATRIX_PARAMETER = "fastCompoundMatrixParameter";

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return FAST_COMPOUND_MATRIX_PARAMETER;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            List<MatrixParameterInterface> matrices = new ArrayList<MatrixParameterInterface>();

            for (int i = 0; i < xo.getChildCount(); ++i) {
                matrices.add((MatrixParameterInterface) xo.getChild(i));
            }

            final String name = xo.hasId() ? xo.getId() : null;

            return new FastCompoundMatrixParameter(name, matrices);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "A fast compound matrix parameter constructed from its component parameters.";
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                new ElementRule(MatrixParameterInterface.class, 1, Integer.MAX_VALUE),
        };

        public Class getReturnType() {
            return FastCompoundMatrixParameter.class;
        }
    };

}

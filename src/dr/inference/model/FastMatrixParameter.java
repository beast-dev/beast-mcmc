/*
 * FastMatrixParameter.java
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

import java.util.ArrayList;
import java.util.List;

/**
 * @author Marc A. Suchard
 */
public class FastMatrixParameter extends CompoundParameter implements MatrixParameterInterface {

    private static final String FAST_MATRIX_PARAMETER = "fastMatrixParameter";
    private static final String ROW_DIMENSION = MatrixParameter.ROW_DIMENSION;
    private static final String COLUMN_DIMENSION = MatrixParameter.COLUMN_DIMENSION;
    private static final String STARTING_VALUE = "startingValue";
    private static final String SIGNAL_COMPONENTS = "signalComponents";

    public FastMatrixParameter(String id, int rowDimension, int colDimension, double startingValue) {
        this(id, rowDimension, colDimension, startingValue, true);
    }

    public FastMatrixParameter(String id, int rowDimension, int colDimension, double startingValue, boolean signalComponents) {
        super(id);
        singleParameter = new Parameter.Default(rowDimension * colDimension);
        addParameter(singleParameter);
        for (int i = 0; i < singleParameter.getDimension(); i++) {
            singleParameter.setParameterValue(i, startingValue);
        }
        this.rowDimension = rowDimension;
        this.colDimension = colDimension;
        DefaultBounds bounds = new DefaultBounds(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, singleParameter.getDimension());
        addBounds(bounds);

        this.signalComponents = signalComponents;
    }

    public FastMatrixParameter(String id, List<Parameter> original, boolean signalComponents) {
        this(id, original.get(0).getDimension(), original.size(), 0.0, signalComponents);

        checkParameterLengths(original);

        setProxyParameterNames(original);


        for (int row = 0; row < rowDimension; ++row) {
            for (int col = 0; col < original.size(); ++col) {
                setParameterValueQuietly(row, col, original.get(col).getParameterValue(row));
            }
        }
        fireParameterChangedEvent(-1, ChangeType.ALL_VALUES_CHANGED);
    }

    private void checkParameterLengths(List<Parameter> parameters) {
        final int length = parameters.get(0).getDimension();
        for (Parameter p : parameters) {
            if (p.getDimension() != length) {
                throw new RuntimeException("All columns must be the same length");
            }
        }
    }

    private void setProxyParameterNames(List<Parameter> original) {
        proxyParameterNames = new ArrayList<String>(original.size());

        for (Parameter p : original) {
            proxyParameterNames.add(p.getParameterName());
        }
    }

    private List<String> proxyParameterNames;

    private String getProxyParameterName(int column) {
        return proxyParameterNames != null ? proxyParameterNames.get(column) : null;
    }

    public Parameter getParameter(int index) {
        if (proxyList == null) {
            proxyList = new ArrayList<ParameterProxy>(colDimension);
            for (int i = 0; i < colDimension; ++i) {
                proxyList.add(new ParameterProxy(this, i));
            }
        }
        return proxyList.get(index);
    }

    public void fireParameterChangedEvent() {
        if (signalComponents) {
            super.fireParameterChangedEvent();
        } else {
            fireParameterChangedEvent(-1, ChangeType.ALL_VALUES_CHANGED);
        }
    }

    class ParameterProxy extends Parameter.Abstract {

        private final int column;
        private final FastMatrixParameter matrix;

        ParameterProxy(FastMatrixParameter matrix, int column) {
            this.matrix = matrix;
//            this.addParameterListener(this.matrix);
            this.column = column;
        }

        @Override
        protected void storeValues() {
            // Do nothing; storeValues() on whole matrix should have been called.
        }

        @Override
        protected void restoreValues() {
            // Do nothing; restoreValues() on whole matrix should have been called.
        }

        @Override
        protected void acceptValues() {
            // Do nothing; acceptValues() on whole matrix should have been called.
        }

        @Override
        protected void adoptValues(Parameter source) {
            throw new RuntimeException("Do not call");
        }

        @Override
        public double getParameterValue(int dim) {
            return matrix.getParameterValue(dim, column);
        }

        @Override
        public void setParameterValue(int dim, double value) {
            matrix.setParameterValue(dim, column, value);
            super.fireParameterChangedEvent(dim, Parameter.ChangeType.VALUE_CHANGED);
        }

        @Override
        public void setParameterValueQuietly(int dim, double value) {
            matrix.setParameterValueQuietly(dim, column, value);
        }

        @Override
        public void setParameterValueNotifyChangedAll(int dim, double value) {
            throw new RuntimeException("Do not call");
        }

        @Override
        public String getParameterName() {
            String proxyName = matrix.getProxyParameterName(column);
            return proxyName != null ? proxyName : getId();
        }

        @Override
        public void addBounds(Bounds<Double> boundary) {
            bounds = boundary;
        }

        @Override
        public Bounds<Double> getBounds() {
            return bounds;
        }

        @Override
        public void fireParameterChangedEvent(int index, ChangeType type) {
            matrix.fireParameterChangedEvent(index + column * getDimension(), type);
            super.fireParameterChangedEvent(index, ChangeType.VALUE_CHANGED);
        }

        @Override
        public void addDimension(int index, double value) {
            throw new RuntimeException("Do not call");
        }

        @Override
        public double removeDimension(int index) {
            throw new RuntimeException("Do not call");
        }

        @Override
        public int getDimension() {
            return matrix.getRowDimension();
        }

        private Bounds<Double> bounds = null;

    }


    @Override
    public double getParameterValue(int row, int col) {
        return singleParameter.getParameterValue(index(row, col));
    }

    @Override
    public double[] getParameterValues() {
        double[] destination = new double[getDimension()];
        copyParameterValues(destination, 0);
        return destination;
    }

    @Override
    public void copyParameterValues(double[] destination, int offset) {
        final double[] source = ((Parameter.Default) singleParameter).inspectParameterValues();
        System.arraycopy(source, 0, destination, offset, source.length);
    }

    @Override
    public void setAllParameterValuesQuietly(double[] values, int offset) {
        final double[] destination = ((Parameter.Default) singleParameter).inspectParameterValues();
        System.arraycopy(values, offset, destination, 0, destination.length);
    }

    @Override
    public String toSymmetricString() {
        return MatrixParameter.toSymmetricString(this);
    }

    @Override
    public boolean isConstrainedSymmetric() {
        return false;
    }

    @Override
    public void setParameterValue(int row, int col, double value) {
        singleParameter.setParameterValue(index(row, col), value);
    }

    @Override
    public void setParameterValueQuietly(int row, int col, double value) {
        singleParameter.setParameterValueQuietly(index(row, col), value);
    }

    @Override
    public void setParameterValueNotifyChangedAll(int row, int col, double value) {
        singleParameter.setParameterValueNotifyChangedAll(index(row, col), value);
    }

    @Override
    public double[] getColumnValues(int col) {
        double[] rtn = new double[getRowDimension()];
        for (int i = 0; i < getRowDimension(); ++i) {
            rtn[i] = getParameterValue(i, col);
        }
        return rtn;
    }

    @Override
    public double[][] getParameterAsMatrix() {
        return MatrixParameterInterface.getParameterAsMatrix(this);
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
    public int getParameterCount() {
        return getColumnDimension();
    }

    @Override
    public int getUniqueParameterCount() {
        return 1;
    }

    @Override
    public Parameter getUniqueParameter(int index) {
        return super.getParameter(0);
    }

    public void addBounds(Bounds<Double> boundary) {
        singleParameter.addBounds(boundary);
    }

    public Bounds<Double> getBounds() {
        return singleParameter.getBounds();
    }


    private final int rowDimension;
    private final int colDimension;
    private final Parameter singleParameter;
    private final boolean signalComponents;

    private List<ParameterProxy> proxyList = null;

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return FAST_MATRIX_PARAMETER;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            final String name = xo.hasId() ? xo.getId() : null;
            final int rowDimension = xo.getIntegerAttribute(ROW_DIMENSION);
            final int colDimension = xo.getIntegerAttribute(COLUMN_DIMENSION);
            final boolean signalComponents = xo.getAttribute(SIGNAL_COMPONENTS, false);

            List<Parameter> columns = new ArrayList<Parameter>();
            for (int i = 0; i < xo.getChildCount(); ++i) {
                if (xo.getChild(i) instanceof Parameter) {
                    columns.add((Parameter) xo.getChild(i));
                }
            }

            if (columns.size() > 0) {
                if (columns.get(0).getDimension() != rowDimension || columns.size() != colDimension) {
                    throw new XMLParseException("Unable to cast matrixParameter to fastMatrixParameter");
                }

                FastMatrixParameter matrix = new FastMatrixParameter(name, columns, signalComponents);
                replaceParameterReferences(xo, matrix);

                return matrix;
            } else {

                final double startingValue = xo.getAttribute(STARTING_VALUE, 1.0);
                return new FastMatrixParameter(name, rowDimension, colDimension, startingValue, signalComponents);

            }
        }

        private void replaceParameterReferences(XMLObject xo, FastMatrixParameter matrix) throws XMLParseException {
            List<XMLObject> children = new ArrayList<XMLObject>();
            for (int i = 0; i < xo.getChildCount(); ++i) {
                if (xo.getChild(i) instanceof Parameter) {
                    Object obj = xo.getRawChild(i);
                    if (obj instanceof Reference) {
                        throw new XMLParseException("Currently can not use references in constructing a fast matrix parameter");
                    }
                    children.add((XMLObject) xo.getRawChild(i));
                }
            }

            if (children.size() != matrix.getColumnDimension()) {
                throw new XMLParseException("Invalid dimensions");
            }

            for (int i = 0; i < children.size(); ++i) {
                children.get(i).setNativeObject(matrix.getParameter(i));
            }
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "A fast matrix parameter constructed from a single parameter.";
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                new ElementRule(Parameter.class, 0, Integer.MAX_VALUE),
                AttributeRule.newIntegerRule(ROW_DIMENSION, false),
                AttributeRule.newIntegerRule(COLUMN_DIMENSION, false),
                AttributeRule.newDoubleRule(STARTING_VALUE, true),
                AttributeRule.newBooleanRule(SIGNAL_COMPONENTS, true),
//                new ElementRule(MatrixParameterInterface.class, true),
        };

        public Class getReturnType() {
            return FastMatrixParameter.class;
        }
    };
}

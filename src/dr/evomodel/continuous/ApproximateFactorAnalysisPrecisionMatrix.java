/*
 * ApproximateFactorAnalysisPrecisionMatrix.java
 *
 * Copyright (c) 2002-2016 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.continuous;

import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import dr.inference.model.*;
import dr.math.matrixAlgebra.Matrix;
import dr.math.matrixAlgebra.RobustSingularValueDecomposition;
import dr.math.matrixAlgebra.Vector;
import dr.xml.*;

/**
 * @author Marc A. Suchard
 * @author Max R. Tolkoff
 *
 */

public class ApproximateFactorAnalysisPrecisionMatrix extends Parameter.Abstract implements
        MatrixParameterInterface, VariableListener {

    private double[] values;
    private double[] storedValues;

    private boolean valuesKnown;
    private boolean storedValuesKnow;

    private final MatrixParameterInterface L;
    private final Parameter gamma;
    private int dim;

    private final static boolean DEBUG = false;

    public ApproximateFactorAnalysisPrecisionMatrix(String name,
                                                    MatrixParameterInterface L,
                                                    Parameter gamma) {
        super(name);

        this.L = L;
        this.gamma = gamma;

        L.addVariableListener(this);
        gamma.addVariableListener(this);

        this.dim = L.getRowDimension();


        values = new double[dim * dim];
    }

    private void computeValues() {
        if (!valuesKnown) {
            computeValuesImp();
            valuesKnown = true;
        }
    }

    private void computeValuesImp() {
        dim = L.getRowDimension();
        double[][] matrix = new double[dim][dim];

        for (int row = 0; row < L.getRowDimension(); ++row) {
            for (int col = 0; col < L.getRowDimension(); ++col) {
                double sum = 0;
                for (int k = 0; k < L.getColumnDimension(); ++k) {
                    sum += L.getParameterValue(row, k) * L.getParameterValue(col, k);
                }
                matrix[row][col] = sum;
            }


        }
//        DoubleMatrix2D LL = new DenseDoubleMatrix2D(matrix);
//        RobustSingularValueDecomposition SVD = new RobustSingularValueDecomposition(LL);
//        double[][] U = SVD.getU().toArray();
//        double[][] V = SVD.getV().toArray();
//        DoubleMatrix2D EvalsTemp= SVD.getS();
//
//        double[] EVals = new double[EvalsTemp.rows()];
//        for (int i = 0; i < EVals.length ; i++) {
//            EVals[i] = EvalsTemp.get(i, i);
//        }
//
//        for (int i = 0; i < EVals.length ; i++) {
//            if(Math.abs(EVals[i]) >= Math.pow(10, -10)){
//                EVals[i] = 1 / EVals[i];
//            }
//            else
//                EVals[i] = 0;
//        }
//
//
//        for (int i = 0; i <U.length; i++) {
//            for (int j = 0; j <V.length ; j++) {
//                matrix[i][j] = 0;
//                for (int k = 0; k < U.length; k++) {
//                    matrix[i][j] += V[i][k] * EVals[k] * U[j][k];
//                }
//
//            }
//
//        }

        for (int row = 0; row < dim; row++) {

            matrix[row][row] += 1 / gamma.getParameterValue(row);
        }

        if (DEBUG) {
            System.err.println("mult:");
            System.err.println(new Matrix(L.getParameterAsMatrix()));
            System.err.println(new Vector(gamma.getParameterValues()) + "\n");
            System.err.println(new Matrix(matrix));
        }

        matrix = new Matrix(matrix).inverse().toComponents();


        int index = 0;
        for (int row = 0; row < dim; ++row) {
            for (int col = 0; col < dim; ++col) {
                values[index] = matrix[row][col];
                ++index;
            }
        }
    }

    @Override
    public String getDimensionName(int index) {
        int row = index % dim + 1;
        int col = index / dim + 1; // column-major
        return getParameterName() + row + col;
    }

    @Override
    public int getDimension() {
        return dim * dim;
    }

    @Override
    public double getParameterValue(int index) {
        computeValues();
        return values[index];
    }

    @Override
    public double[][] getParameterAsMatrix() {
        computeValues();
        double[][] matrix = new double[dim][dim];
        for (int i = 0; i < dim; ++i) {
            System.arraycopy(values, i * dim, matrix[i], 0, dim);
        }

        if (DEBUG) {
            System.err.println("vec:");
            System.err.println(new Vector(values));
            System.err.println(new Matrix(matrix));
         System.err.println("");
        }

        return matrix;
    }

    @Override
    public double getParameterValue(int row, int col) {
        computeValues();
        return values[col * dim + row]; // column-major
    }

    @Override
    public Parameter getParameter(int column) {
        return null;
    }

    @Override
    public double[] getParameterValues() {
        computeValues();
        double[] matrix = new double[values.length];
        System.arraycopy(values, 0, matrix, 0, values.length);
        return matrix;
    }

    @Override
    protected void storeValues() {
        L.storeParameterValues();
        gamma.storeParameterValues();

        if (storedValues == null) {
            storedValues = new double[dim * dim];
        }

        System.arraycopy(values, 0, storedValues, 0, values.length);
        storedValuesKnow = valuesKnown;
    }

    @Override
    protected void restoreValues() {
        L.restoreParameterValues();
        gamma.restoreParameterValues();

        double[] tmp = values;
        values = storedValues;
        storedValues = tmp;
        valuesKnown = storedValuesKnow;
    }

    @Override
    protected void acceptValues() {
        L.acceptParameterValues();
        gamma.acceptParameterValues();
    }

    @Override
    protected void adoptValues(Parameter source) {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public void setParameterValue(int dim, double value) {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public void setParameterValueQuietly(int dim, double value) {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public void setParameterValueNotifyChangedAll(int dim, double value) {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public String getParameterName() {
        return getId();
    }

    @Override
    public void addBounds(Bounds<Double> bounds) {

    }

    @Override
    public Bounds<Double> getBounds() {
        return null;
    }

    @Override
    public void addDimension(int index, double value) {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public double removeDimension(int index) {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public void setParameterValue(int row, int col, double value) {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public void setParameterValueQuietly(int row, int col, double value) {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public void setParameterValueNotifyChangedAll(int row, int col, double value) {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public double[] getColumnValues(int col) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public int getColumnDimension() {
        return dim;
    }

    @Override
    public int getRowDimension() {
        return dim;
    }

    @Override
    public int getUniqueParameterCount() {
        return 2;
    }

    @Override
    public Parameter getUniqueParameter(int index) {
        return (index == 0) ? L : gamma;
    }

    @Override
    public void copyParameterValues(double[] destination, int offset) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public void setAllParameterValuesQuietly(double[] values, int offset) {
        throw new RuntimeException("Not implemented");
    }

    @Override
    public String toSymmetricString() {
        return MatrixParameter.toSymmetricString(this);
    }

    @Override
    public void variableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        valuesKnown = false;
        fireParameterChangedEvent();
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public static final String APPROXIMATE_PARAMETER = "approximateFactorAnalysisPrecision";
        public static final String L_LABEL = "L";
        public static final String GAMMA = "gamma";

        public String getParserName() {
            return APPROXIMATE_PARAMETER;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            MatrixParameterInterface L = (MatrixParameterInterface) xo.getElementFirstChild(L_LABEL);
            Parameter gamma = (Parameter) xo.getElementFirstChild(GAMMA);

            String name = xo.hasId() ? xo.getId() : APPROXIMATE_PARAMETER;

            return new ApproximateFactorAnalysisPrecisionMatrix(name, L, gamma);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "A diffusion approximation to a factor analysis";
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                new ElementRule(L_LABEL, new XMLSyntaxRule[] {
                        new ElementRule(MatrixParameterInterface.class),
                }),
                new ElementRule(GAMMA, new XMLSyntaxRule[] {
                        new ElementRule(Parameter.class),
                }),
        };

        public Class getReturnType() {
            return ApproximateFactorAnalysisPrecisionMatrix.class;
        }
    };
}
//        }
//                DoubleMatrix2D LL = new DenseDoubleMatrix2D(matrix);
//                RobustSingularValueDecomposition SVD = new RobustSingularValueDecomposition(LL);
//                double[][] U = SVD.getU().toArray();
//                double[][] V = SVD.getV().toArray();
//                double[] EVals = SVD.getSingularValues();
//
//                for (int i = 0; i <EVals.length ; i++) {
//        if(EVals[i] != 0){
//        EVals[i] = 1 / EVals[i];
//        }
//        }
//        for (int i = 0; i < EVals.length; i++) {
//        for (int j = 0; j < EVals.length; j++) {
//        U[i][j] = U[i][j] * EVals[i];
//        }
//        }
//
//        for (int i = 0; i <U.length; i++) {
//        for (int j = 0; j <V.length ; j++) {
//        matrix[i][j] = U[j][i] * V[i][j];
//        }
//
//        }
//
//        for (int row = 0; row < dim; row++) {
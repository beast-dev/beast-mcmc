/*
 * MarkovRandomFieldMatrix.java
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

import dr.util.Transform;

/**
 * @author Marc Suchard
 */
public class MarkovRandomFieldMatrix extends MatrixParameter {

    private Parameter diagonalParameter;
    private Parameter offDiagonalParameter;

    private final Parameter nuggetParameter;

    private final Transform diagonalTransform;
    private final Transform offDiagonalTransform;

    private boolean asCorrelation = false;

    private int dim;

//    public MarkovRandomFieldMatrix(String name, int dim, Parameter diagonals, Parameter offDiagonal,
//                                   Parameter nugget,
//                                   boolean asCorrelation) {
//        this(name, dim, diagonals, offDiagonal, nugget, asCorrelation, null, null);
//    }

    public MarkovRandomFieldMatrix(String name, int dim, Parameter diagonals, Parameter offDiagonal,
                                   Parameter nugget,
                                   boolean asCorrelation,
                                   Transform diagonalTransform, Transform offDiagonalTransform) {
        super(name);
        diagonalParameter = diagonals;
        offDiagonalParameter = offDiagonal;
        nuggetParameter = nugget;

        addParameter(diagonalParameter);
        addParameter(offDiagonalParameter);
        addParameter(nuggetParameter);

        this.dim = dim;
        this.asCorrelation = asCorrelation;

        this.diagonalTransform = (diagonalTransform != null) ? diagonalTransform : Transform.NONE;
        this.offDiagonalTransform = (offDiagonalTransform != null) ? offDiagonalTransform : Transform.NONE;
    }

//    public void variableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
//        if (variable == diagonalParameter || variable == offDiagonalParameter || variable == nuggetParameter) {
//            fireParameterChangedEvent(-1, ChangeType.ALL_VALUES_CHANGED);
//        } else {
//            throw new IllegalArgumentException("Unknown variable '" + variable.getVariableName() + "'");
//        }
//    }

//    public double[] getAttributeValue() {
//        double[] stats = new double[dim * dim];
//        int index = 0;
//        for (int i = 0; i < dim; i++) {
//            for (int j = 0; j < dim; j++) {
//                stats[index] = getParameterValue(i, j);
//                index++;
//            }
//        }
//        return stats;
//    }
//
    public int getDimension() {
        return dim * dim;
    }

    public String getDimensionName(int i) {
        int row = i / dim;
        int col = i % dim;
        return getParameterName() + "_" + (row + 1) + "_" + (col + 1);
    }

    private double getDiagonalParameterValue(int i) {
        return diagonalTransform.inverse(diagonalParameter.getParameterValue(i));
    }

    private double getOffDiagonalParameterValue(int i) {
        return offDiagonalTransform.inverse(offDiagonalParameter.getParameterValue(i));
    }

    private double getNuggetValue(int i) {
        return nuggetParameter.getParameterValue(i);
    }

    public double getParameterValue(int i) {
        int row = i / dim;
        int col = i % dim;
        return getParameterValue(row, col);
    }

    public double getParameterValue(int row, int col) {
        if (row == col) {
            double diag = getDiagonalParameterValue(0);
            if (row > 0 && row < (dim - 1)) {
                diag *= 2; // TODO Assumes RW-1 model
            }
            diag += getNuggetValue(0);
            return diag;
//            return getDiagonalParameterValue(row);
        } else if (row == (col - 1) || row == (col + 1)) {
            if (asCorrelation) {
                return -getOffDiagonalParameterValue(0) *
                        Math.sqrt(getDiagonalParameterValue(0) * getDiagonalParameterValue(0));
            }
            return getOffDiagonalParameterValue(0);
        }
        return 0.0;
    }

//    public double[][] getParameterAsMatrix() {
//        final int I = dim;
//        double[][] parameterAsMatrix = new double[I][I];
//        for (int i = 0; i < I; i++) {
//            parameterAsMatrix[i][i] = getParameterValue(i, i);
//            for (int j = i + 1; j < I; j++) {
//                parameterAsMatrix[j][i] = parameterAsMatrix[i][j] = getParameterValue(i, j);
//            }
//        }
//        return parameterAsMatrix;
//    }

    public int getColumnDimension() {
        return dim;
    }

    public int getRowDimension() {
        return dim;
    }
}

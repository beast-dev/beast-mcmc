/*
 * WishartSufficientStatistics.java
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

package dr.math.distributions;

import dr.math.matrixAlgebra.Vector;

import java.util.Arrays;

public class WishartSufficientStatistics {

    public WishartSufficientStatistics(int dim) {
        df = 0;
        scaleMatrix = new double[dim * dim];
    }

    public WishartSufficientStatistics(int df, double[] scaleMatrix) {
        this.df = df;
        this.scaleMatrix = scaleMatrix;
    }

    public WishartSufficientStatistics clone() {
        return new WishartSufficientStatistics(this.df, scaleMatrix.clone());
    }

    public void copyTo(WishartSufficientStatistics destination) {
        assert (destination != null);

        destination.df = this.df;
        System.arraycopy(this.scaleMatrix, 0, destination.scaleMatrix, 0, this.scaleMatrix.length);
    }

//    public WishartSufficientStatistics(int df, double[] matrix, int dim) {
//        this(df, convertToSquare(matrix, dim));
//    }
//
//    public WishartSufficientStatistics(int df, double[] matrix, int dim, int count) {
//        this(df, convertToSquareAndAccumulate(matrix, dim, count));
//    }
//
//    public WishartSufficientStatistics(double[] matrix, int dim, int count) {
//        this(convertFromSquareAndAccumulate(matrix, dim, count));
//    }

    public WishartSufficientStatistics(int[] dfs, double[] matrix) {
        this(sum(dfs), accumulateMatrix(matrix, dfs.length));
    }

//    private WishartSufficientStatistics(WishartSufficientStatistics copy) {
//        this.df = copy.df;
//        this.scaleMatrix = copy.scaleMatrix;
//    }

    public final int getDf() {
        return df;
    }

    public final double[] getScaleMatrix() {
        return scaleMatrix;
    }
    
    public final void incrementDf(int n) {
        df += n;
    }

    public final void clear() {
        df = 0;
        Arrays.fill(scaleMatrix, 0.0);
    }

//    private static double[][] convertToSquare(double[] vector, int dim)  {
//        double[][] matrix = new double[dim][dim];
//        int offset = 0;
//        for (int row = 0; row < dim; ++row) {
//            System.arraycopy(vector, offset, matrix[row], 0, dim);
//            offset += dim;
//        }
//        return matrix;
//    }

//    private static WishartSufficientStatistics convertFromSquareAndAccumulate(double[] vector, int dim, int count)  {
//        double[][] matrix = new double[dim][dim];
//        int df = 0;
//
//        for (int c = 0; c < count; ++c) {
//            int offset = (dim * dim + 1) * c;
//            for (int i = 0; i < dim; ++i) {
//                for (int j = 0; j < dim; ++j) {
//                    matrix[i][j] += vector[offset];
//                    ++offset;
//                }
//            }
//            df += vector[offset];
//        }
//        return new WishartSufficientStatistics(df, matrix);
//    }

    private static double[] accumulateMatrix(double[] matrix, final int count) {
        if (count == 1) {
            return matrix;
        } else {
            final int length = matrix.length / count;
            double[] result = new double[length];
            System.arraycopy(matrix, 0, result, 0, length);

            for (int i = 1; i < count; ++i) {
                final int offset = i * length;
                for (int j = 0; j < length; ++j) {
                    result[j] += matrix[offset + j];
                }
            }
            return result;
        }
    }

    private static int sum(int[] vector) {
        int sum = 0;
        for (int x : vector) {
            sum += x;
        }
        return sum;
    }

//    private static double[][] convertToSquareAndAccumulate(double[] vector, int dim, int count)  {
//        double[][] matrix = new double[dim][dim];
//
//        for (int c = 0; c < count; ++c) {
//            int offset = dim * dim * c;
//            for (int i = 0; i < dim; ++i) {
//                for (int j = 0; j < dim; ++j) {
//                    matrix[i][j] += vector[offset];
//                    ++offset;
//                }
//            }
//        }
//        return matrix;
//    }

    public String toString() {
        return df + " : " + new Vector(scaleMatrix).toString();
    }

    private int df;
    private final double[] scaleMatrix;
}
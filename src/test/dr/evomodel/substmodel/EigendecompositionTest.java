/*
 * EigendecompositionTest.java
 *
 * Copyright (c) 2002-2012 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package test.dr.evomodel.substmodel;

import dr.app.beagle.evomodel.substmodel.ComplexSubstitutionModel;
import dr.app.beagle.evomodel.substmodel.EigenDecomposition;
import dr.app.beagle.evomodel.substmodel.FrequencyModel;
import dr.evolution.datatype.GeneralDataType;
import dr.inference.model.Parameter;
import dr.math.MathUtils;
import test.dr.math.MathTestCase;

import java.util.ArrayList;
import java.util.Arrays;

public class EigendecompositionTest extends MathTestCase {

    int dim;
    int[][] columnNumberLookup;
    double tolerance;

    public EigendecompositionTest(int dim, double tolerance) {
        this.dim = dim;
        this.tolerance = tolerance;
        columnNumberLookup = new int[dim][dim];
        int columnCount = 0;
        for (int row = 0; row < dim; row++) {
            for (int column = 0; column < dim; column++) {
                if (row == column) {
                    columnNumberLookup[row][column] = -1;
                } else if (row < column) {
                    columnNumberLookup[row][column] = columnCount;
                    columnCount++;
                }
            }
        }
        for (int column = 0; column < dim; column++) {
            for (int row = 0; row < dim; row++) {
                if (row > column) {
                    columnNumberLookup[row][column] = columnCount;
                    columnCount++;
                }
            }
        }
    }

    public void doTest() {

        String[] labels = new String[dim];
        for (Integer i = 0; i < dim; i++) {
            labels[i] = i.toString();
        }
        GeneralDataType dataType = new GeneralDataType(labels);
        Parameter.Default freqVector = new Parameter.Default(dim);
        for (int i = 0; i < dim; i++) {
            freqVector.setParameterValue(i, (double) 1 / dim);
        }
        FrequencyModel freqModel = new FrequencyModel(dataType, freqVector);
        Parameter.Default matrixParameter = new Parameter.Default(dim * (dim - 1));
        ComplexSubstitutionModel substModel = new ComplexSubstitutionModel("test", dataType, freqModel, matrixParameter);

        ArrayList<Integer[]> indicators = getAllPossibleZeroRows(dim, columnNumberLookup);

        for (Integer[] ind : indicators) {

            int successED = 0;
            int[] successTP = new int[dim];

            int[] zeroRows = new int[dim];


            for (int row = 0; row < dim; row++) {
                boolean zeroRow = true;
                for (int col = 0; col < dim; col++) {
                    if (row != col && ind[columnNumberLookup[row][col]] != 0) {
                        zeroRow = false;
                    }
                }
                zeroRows[row] = zeroRow ? 1 : 0;
            }

            for (int run = 0; run < 100; run++) {

                for (int i = 0; i < dim * (dim - 1); i++) {
                    matrixParameter.setParameterValue(i, MathUtils.nextDouble() * ind[i]);
                }

                normalise(matrixParameter);

                boolean EigenDecompOK = testEigenDecomposition(substModel);

                Boolean[] transProbsOK = new Boolean[dim];

                for (int i = 0; i < dim; i++) {
                    transProbsOK[i] = zeroRows[i] != 1 || testProbMatrix(substModel, i, dim, tolerance);
                }

                if (EigenDecompOK) {
                    successED++;
                }
                for (int i = 0; i < dim; i++)
                    if (EigenDecompOK && transProbsOK[i]) {
                        successTP[i]++;
                    }
            }

            System.out.print("Indicators: ");
            for (int i = 0; i < dim * (dim - 1); i++) {
                System.out.print(ind[i] + " ");
            }
            System.out.println();
            System.out.println("Eigendecomposition succeeded " + successED + " times out of 100");
            for (int i = 0; i < dim; i++) {
                if (zeroRows[i] == 1) {
                    System.out.println("Row " + (i + 1) + ": Connectivity check would work " + successTP[i] + " times out of " + successED);
                }
            }
            System.out.println();
        }

    }

    private static boolean testEigenDecomposition(ComplexSubstitutionModel model) {
        boolean out = false;
        EigenDecomposition eigen = model.getEigenDecomposition();

        for (int i = 0; i < eigen.getEigenVectors().length; i++) {
            if (eigen.getEigenVectors()[i] != 0) {
                out = true;
            }
        }

        return out;

    }

    private static ArrayList<Integer[]> getAllPossibleZeroRows(int dim, int[][] lookup) {
        ArrayList<Integer[]> allPossibleZeroRows = new ArrayList<Integer[]>();
        for (int zeroRow = 0; zeroRow < dim; zeroRow++) {
            Integer[] possibility = new Integer[dim * (dim - 1)];
            for (int row = 0; row < dim; row++) {
                for (int col = 0; col < dim; col++) {
                    if (row != col) {
                        possibility[lookup[row][col]] = row == zeroRow ? 0 : 1;
                    }
                }
            }
            allPossibleZeroRows.add(possibility);
        }
        return allPossibleZeroRows;
    }

    private static ArrayList<Integer[]> getAllPossibilities(int length) {
        ArrayList<Integer[]> allPossibilities = new ArrayList<Integer[]>();
        if (length == 1) {
            Integer[] one = new Integer[1];
            one[0] = 1;
            Integer[] zero = new Integer[1];
            zero[0] = 0;
            allPossibilities.add(one);
            allPossibilities.add(zero);
        } else {
            for (Integer[] tail : getAllPossibilities(length - 1)) {
                Integer[] one = new Integer[1];
                one[0] = 1;
                Integer[] zero = new Integer[1];
                zero[0] = 0;
                Integer[] onePlus = concat(one, tail);
                Integer[] zeroPlus = concat(zero, tail);
                allPossibilities.add(onePlus);
                allPossibilities.add(zeroPlus);
            }
        }
        return allPossibilities;
    }

    public static <T> T[] concat(T[] first, T[] second) {
        T[] result = Arrays.copyOf(first, first.length + second.length);
        System.arraycopy(second, 0, result, first.length, second.length);
        return result;
    }

    private static Parameter.Default normalise(Parameter.Default matrix) {
        double sum = getMatrixSum(matrix);
        int dim = (int) Math.sqrt(matrix.getDimension());

        for (int i = 0; i < dim * (dim - 1); i++) {
            matrix.setParameterValue(i, matrix.getParameterValue(i) * dim / sum);
        }

        return matrix;
    }

    private static boolean testProbMatrix(ComplexSubstitutionModel model, int rowOfInterest, int dim, double tolerance) {
        boolean out = true;

        double[] probs = new double[dim * dim];

        model.getTransitionProbabilities(1, probs);

        for (int i = rowOfInterest * dim; i < rowOfInterest * dim + dim; i++) {
            if (i != rowOfInterest * dim + rowOfInterest) {
                if (probs[i] > tolerance) {
                    out = false;
                }
            }
        }
        return out;

    }

    public static double getMatrixSum(Parameter.Default matrix) {
        double sum = 0;
        for (int i = 0; i < matrix.getDimension(); i++) {
            sum += matrix.getParameterValue(i);
        }
        return sum;
    }

    public static void main(String[] args) {
        EigendecompositionTest run = new EigendecompositionTest(4, 1E-20);
        run.doTest();
    }

}


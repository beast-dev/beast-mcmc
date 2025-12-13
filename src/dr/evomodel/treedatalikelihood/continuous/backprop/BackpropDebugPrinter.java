package dr.evomodel.treedatalikelihood.continuous.backprop;

/*
 * BackpropDebugPrinter.java
 *
 * Utility class for debug printing in backpropagation gradient computation.
 * Provides formatted output for matrices, vectors, and arrays with smart
 * handling of both small and large matrices.
 *
 *  BEAST is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 */

import org.ejml.data.DenseMatrix64F;

/**
 * Debug printing utilities for backpropagation.
 *
 * Features:
 *   - Smart matrix printing (full for small, corners+stats for large)
 *   - Consistent formatting across all debug output
 *   - Null-safe printing
 *   - Statistics computation (range, norms)
 */
public class BackpropDebugPrinter {

    private static final int SMALL_MATRIX_THRESHOLD = 5;
    private static final int CORNER_SIZE = 3;
    private static final int MAX_ARRAY_DISPLAY = 20;

    // ------------------------------------------------------------------------
    // Matrix printing
    // ------------------------------------------------------------------------

    /**
     * Print a matrix with automatic format selection.
     * Small matrices (≤5×5) are printed in full.
     * Large matrices show corners and statistics.
     */
    public static void printMatrix(String name, DenseMatrix64F mat) {
        if (mat == null) {
            System.err.println(name + ": null");
            return;
        }

        System.err.println(name + " (" + mat.numRows + "×" + mat.numCols + "):");

        if (isSmall(mat)) {
            printFullMatrix(mat);
        } else {
            printLargeMatrix(mat);
        }
    }

    /**
     * Print matrix with custom value formatter.
     */
    public static void printMatrix(String name, DenseMatrix64F mat, String format) {
        if (mat == null) {
            System.err.println(name + ": null");
            return;
        }

        System.err.println(name + " (" + mat.numRows + "×" + mat.numCols + "):");

        if (isSmall(mat)) {
            printFullMatrix(mat, format);
        } else {
            printLargeMatrix(mat);
        }
    }

    /**
     * Print full matrix (for small matrices).
     */
    private static void printFullMatrix(DenseMatrix64F mat) {
        printFullMatrix(mat, "%12.6f");
    }

    private static void printFullMatrix(DenseMatrix64F mat, String format) {
        for (int i = 0; i < mat.numRows; i++) {
            System.err.print("  [");
            for (int j = 0; j < mat.numCols; j++) {
                System.err.printf(format, mat.get(i, j));
                if (j < mat.numCols - 1) System.err.print(", ");
            }
            System.err.println("]");
        }
    }

    /**
     * Print large matrix (corners + statistics).
     */
    private static void printLargeMatrix(DenseMatrix64F mat) {
        System.err.println("  (Large matrix - showing corners only)");

        // Top-left corner
        System.err.println("  Top-left:");
        int rowsToShow = Math.min(CORNER_SIZE, mat.numRows);
        int colsToShow = Math.min(CORNER_SIZE, mat.numCols);

        for (int i = 0; i < rowsToShow; i++) {
            System.err.print("    [");
            for (int j = 0; j < colsToShow; j++) {
                System.err.printf("%12.6f", mat.get(i, j));
                if (j < colsToShow - 1) System.err.print(", ");
            }
            if (mat.numCols > CORNER_SIZE) {
                System.err.print(", ...");
            }
            System.err.println("]");
        }

        if (mat.numRows > CORNER_SIZE) {
            System.err.println("  ...");
        }

        // Statistics
        MatrixStats stats = computeStats(mat);
        System.err.printf("  Range: [%12.6f, %12.6f]\n", stats.min, stats.max);
        System.err.printf("  Frobenius norm: %12.6f\n", stats.frobeniusNorm);
        System.err.printf("  Mean: %12.6f, Std dev: %12.6f\n", stats.mean, stats.stdDev);
    }

    /**
     * Check if matrix is "small" (should be printed in full).
     */
    private static boolean isSmall(DenseMatrix64F mat) {
        return mat.numRows <= SMALL_MATRIX_THRESHOLD && mat.numCols <= SMALL_MATRIX_THRESHOLD;
    }

    // ------------------------------------------------------------------------
    // Array printing
    // ------------------------------------------------------------------------

    /**
     * Print a double array with optional truncation.
     */
    public static void printDoubleArray(String name, double[] arr) {
        if (arr == null) {
            System.err.println(name + ": null");
            return;
        }

        System.err.print(name + " (length " + arr.length + "): [");

        int displayLength = Math.min(arr.length, MAX_ARRAY_DISPLAY);
        for (int i = 0; i < displayLength; i++) {
            System.err.printf("%.6f", arr[i]);
            if (i < displayLength - 1) System.err.print(", ");
        }

        if (arr.length > MAX_ARRAY_DISPLAY) {
            System.err.print(", ... (" + (arr.length - MAX_ARRAY_DISPLAY) + " more)");
        }
        System.err.println("]");
    }

    /**
     * Print array with custom format.
     */
    public static void printDoubleArray(String name, double[] arr, String format) {
        if (arr == null) {
            System.err.println(name + ": null");
            return;
        }

        System.err.print(name + " (length " + arr.length + "): [");

        int displayLength = Math.min(arr.length, MAX_ARRAY_DISPLAY);
        for (int i = 0; i < displayLength; i++) {
            System.err.printf(format, arr[i]);
            if (i < displayLength - 1) System.err.print(", ");
        }

        if (arr.length > MAX_ARRAY_DISPLAY) {
            System.err.print(", ... (" + (arr.length - MAX_ARRAY_DISPLAY) + " more)");
        }
        System.err.println("]");
    }

    // ------------------------------------------------------------------------
    // Statistics computation
    // ------------------------------------------------------------------------

    /**
     * Container for matrix statistics.
     */
    public static class MatrixStats {
        public final double min;
        public final double max;
        public final double mean;
        public final double stdDev;
        public final double frobeniusNorm;

        public MatrixStats(double min, double max, double mean, double stdDev, double frobeniusNorm) {
            this.min = min;
            this.max = max;
            this.mean = mean;
            this.stdDev = stdDev;
            this.frobeniusNorm = frobeniusNorm;
        }
    }

    /**
     * Compute statistics for a matrix.
     */
    public static MatrixStats computeStats(DenseMatrix64F mat) {
        double min = Double.POSITIVE_INFINITY;
        double max = Double.NEGATIVE_INFINITY;
        double sum = 0.0;
        double sumSq = 0.0;
        int count = mat.numRows * mat.numCols;

        for (int i = 0; i < mat.numRows; i++) {
            for (int j = 0; j < mat.numCols; j++) {
                double val = mat.get(i, j);
                if (val < min) min = val;
                if (val > max) max = val;
                sum += val;
                sumSq += val * val;
            }
        }

        double mean = sum / count;
        double variance = (sumSq / count) - (mean * mean);
        double stdDev = Math.sqrt(Math.max(0.0, variance));  // Avoid numerical issues
        double frobeniusNorm = Math.sqrt(sumSq);

        return new MatrixStats(min, max, mean, stdDev, frobeniusNorm);
    }

    /**
     * Compute L2 norm of a vector/array.
     */
    public static double computeL2Norm(double[] arr) {
        double sum = 0.0;
        for (double v : arr) {
            sum += v * v;
        }
        return Math.sqrt(sum);
    }

    /**
     * Compute L2 norm of a matrix (Frobenius norm).
     */
    public static double computeL2Norm(DenseMatrix64F mat) {
        double sum = 0.0;
        for (int i = 0; i < mat.numRows; i++) {
            for (int j = 0; j < mat.numCols; j++) {
                double val = mat.get(i, j);
                sum += val * val;
            }
        }
        return Math.sqrt(sum);
    }

    // ------------------------------------------------------------------------
    // Specialized printing
    // ------------------------------------------------------------------------

    /**
     * Print a section header.
     */
    public static void printSectionHeader(String title) {
        int totalWidth = 60;
        int padding = (totalWidth - title.length() - 2) / 2;

        System.err.println();
        System.err.println("=".repeat(totalWidth));
        System.err.println(" ".repeat(padding) + title);
        System.err.println("=".repeat(totalWidth));
    }

    /**
     * Print a subsection header.
     */
    public static void printSubsectionHeader(String title) {
        System.err.println();
        System.err.println("--- " + title + " ---");
    }

    /**
     * Print a key-value pair.
     */
    public static void printKeyValue(String key, Object value) {
        System.err.printf("%s: %s\n", key, value);
    }

    /**
     * Print a scalar with label.
     */
    public static void printScalar(String name, double value) {
        System.err.printf("%s = %.6f\n", name, value);
    }

    /**
     * Print a scalar with custom format.
     */
    public static void printScalar(String name, double value, String format) {
        System.err.printf("%s = " + format + "\n", name, value);
    }

    // ------------------------------------------------------------------------
    // Comparison utilities
    // ------------------------------------------------------------------------

    /**
     * Print comparison between two matrices.
     */
    public static void printMatrixComparison(String name,
                                             DenseMatrix64F mat1, String label1,
                                             DenseMatrix64F mat2, String label2) {
        System.err.println("\n--- Comparison: " + name + " ---");

        if (mat1 == null || mat2 == null) {
            System.err.println("One or both matrices are null");
            return;
        }

        if (mat1.numRows != mat2.numRows || mat1.numCols != mat2.numCols) {
            System.err.println("Dimension mismatch!");
            System.err.println(label1 + ": " + mat1.numRows + "×" + mat1.numCols);
            System.err.println(label2 + ": " + mat2.numRows + "×" + mat2.numCols);
            return;
        }

        // Compute difference
        DenseMatrix64F diff = new DenseMatrix64F(mat1.numRows, mat1.numCols);
        for (int i = 0; i < mat1.numRows; i++) {
            for (int j = 0; j < mat1.numCols; j++) {
                diff.set(i, j, mat1.get(i, j) - mat2.get(i, j));
            }
        }

        MatrixStats stats1 = computeStats(mat1);
        MatrixStats stats2 = computeStats(mat2);
        MatrixStats statsDiff = computeStats(diff);

        System.err.printf("%s - Range: [%.6f, %.6f], Norm: %.6f\n",
                label1, stats1.min, stats1.max, stats1.frobeniusNorm);
        System.err.printf("%s - Range: [%.6f, %.6f], Norm: %.6f\n",
                label2, stats2.min, stats2.max, stats2.frobeniusNorm);
        System.err.printf("Difference - Max abs: %.6e, RMS: %.6e\n",
                Math.max(Math.abs(statsDiff.min), Math.abs(statsDiff.max)),
                statsDiff.frobeniusNorm / Math.sqrt(diff.numRows * diff.numCols));
    }

    /**
     * Print comparison between two arrays.
     */
    public static void printArrayComparison(String name,
                                            double[] arr1, String label1,
                                            double[] arr2, String label2) {
        System.err.println("\n--- Comparison: " + name + " ---");

        if (arr1 == null || arr2 == null) {
            System.err.println("One or both arrays are null");
            return;
        }

        if (arr1.length != arr2.length) {
            System.err.println("Length mismatch!");
            System.err.println(label1 + ": length " + arr1.length);
            System.err.println(label2 + ": length " + arr2.length);
            return;
        }

        double maxDiff = 0.0;
        double sumSqDiff = 0.0;

        for (int i = 0; i < arr1.length; i++) {
            double diff = Math.abs(arr1[i] - arr2[i]);
            maxDiff = Math.max(maxDiff, diff);
            sumSqDiff += diff * diff;
        }

        double rmsDiff = Math.sqrt(sumSqDiff / arr1.length);

        System.err.printf("%s norm: %.6f\n", label1, computeL2Norm(arr1));
        System.err.printf("%s norm: %.6f\n", label2, computeL2Norm(arr2));
        System.err.printf("Max abs difference: %.6e\n", maxDiff);
        System.err.printf("RMS difference: %.6e\n", rmsDiff);
    }
}
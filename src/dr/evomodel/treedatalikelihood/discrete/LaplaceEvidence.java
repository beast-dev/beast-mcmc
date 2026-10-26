/*
 * LaplaceEvidence.java
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

package dr.evomodel.treedatalikelihood.discrete;

import com.github.lbfgs4j.liblbfgs.Function;

import java.util.ArrayList;
import java.util.List;

/**
 * The curvature of a maximizer's objective at its optimum, in the coordinates it optimizes (y, after the transform),
 * and the Laplace approximation of the log evidence (plan 8 X3), shared by LaplaceApproximation and
 * ApproximateTreeDataLikelihood.
 * <p>
 * With the transformed density l(y) = log p(x(y)) - logJacobian(x(y)) (BEAST's convention: logJacobian is log|dy/dx|,
 * the density HMC samples), maximized at y, and A = -d2 l/dy2 there,
 * <pre>
 *     log Z = l(y) + (k/2) log 2 pi - (1/2) log det A,
 * </pre>
 * which for a Gaussian of standard deviation s gives (1/2) log 2 pi + log s. It is the mode of the transformed density
 * that the approximation needs, so the evidence requires a maximizer that includes the Jacobian; the curvature
 * (and the covariance, for Wald intervals in y) does not.
 * <p>
 * Coordinates held at a bound and the stick-breaking simplex's dependent coordinate are left out: the evidence and
 * covariance are conditional on the held values. The Hessian is the central difference of the objective's gradient,
 * symmetrized; its log determinant comes from a Cholesky factorization, so a matrix that is not positive definite
 * (not a maximum, or flat in some direction) fails instead of giving a number.
 *
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public class LaplaceEvidence {

    private static final double RELATIVE_STEP = 1E-4;

    private final MaximizerWrtParameter maximizer;
    private final int[] coordinates; // the free coordinates of the optimum
    private final double[] mode; // the optimum, in y
    private final double[][] precision; // A, over the free coordinates
    private final double logDet; // log det A
    private final double[][] covariance; // A^-1

    /**
     * Computes the curvature at the maximizer's latest optimum.
     *
     * @throws IllegalStateException when the maximizer has not run, or A is not positive definite
     */
    public LaplaceEvidence(MaximizerWrtParameter maximizer) {
        if (!maximizer.wasExecuted()) {
            throw new IllegalStateException("The maximizer has not run: there is no optimum to approximate at");
        }
        this.maximizer = maximizer;
        this.mode = maximizer.getMinimumPoint(false).clone();

        List<Integer> free = new ArrayList<Integer>();
        List<MaximizerWrtParameter.Coordinate> optimum = maximizer.getOptimum();
        for (int i = 0; i < optimum.size(); ++i) {
            if (optimum.get(i).free && !optimum.get(i).held) {
                free.add(i);
            }
        }
        coordinates = new int[free.size()];
        for (int i = 0; i < coordinates.length; ++i) {
            coordinates[i] = free.get(i);
        }

        precision = negativeHessian(maximizer.getFunction(), mode, coordinates);
        double[][] lower = cholesky(precision);
        double sum = 0.0;
        for (int i = 0; i < lower.length; ++i) {
            sum += Math.log(lower[i][i]);
        }
        logDet = 2.0 * sum;
        covariance = inverse(lower);
    }

    /**
     * @return log Z, the Laplace approximation of the log evidence: defined only when the maximizer maximized the
     * transformed density (includeJacobian), whose mode the approximation is taken at
     * @throws IllegalStateException otherwise
     */
    public double getLogEvidence() {
        if (!maximizer.includesJacobian() && maximizer.getTransform() != null) {
            throw new IllegalStateException("The Laplace evidence needs the mode of the transformed density: maximize " +
                    "with includeJacobian=\"true\"");
        }
        // the objective F is -l(y)
        return -maximizer.getMinimumValue() + coordinates.length / 2.0 * Math.log(2.0 * Math.PI) - 0.5 * logDet;
    }

    /**
     * @return the indices, in the maximized parameter, of the coordinates the curvature covers
     */
    public int[] getCoordinates() {
        return coordinates.clone();
    }

    /**
     * @return the optimum in y, over all coordinates
     */
    public double[] getMode() {
        return mode.clone();
    }

    /**
     * @return A, the negative Hessian of the objective's log density in y, over getCoordinates()
     */
    public double[][] getPrecision() {
        return copy(precision);
    }

    public double[][] getCovariance() {
        return copy(covariance);
    }

    public double getLogDeterminant() {
        return logDet;
    }

    // ---------------------------------------------------------------------------------------------------------------

    // -d2 l/dy2 over the coordinates, by central differences of the objective's gradient (that of -l), symmetrized;
    // the model is left at the mode
    static double[][] negativeHessian(Function function, double[] mode, int[] coordinates) {
        final int k = coordinates.length;
        final double[][] hessian = new double[k][k];
        final double[] point = mode.clone();
        for (int b = 0; b < k; ++b) {
            final int j = coordinates[b];
            final double h = RELATIVE_STEP * Math.max(1.0, Math.abs(mode[j]));
            point[j] = mode[j] + h;
            final double[] plus = function.gradientAt(point);
            point[j] = mode[j] - h;
            final double[] minus = function.gradientAt(point);
            point[j] = mode[j];
            for (int a = 0; a < k; ++a) {
                final int i = coordinates[a];
                // the gradient of the objective is -dl/dy, so this is -d2 l/dy_i dy_j
                hessian[a][b] = (plus[i] - minus[i]) / (2.0 * h);
            }
        }
        function.valueAt(mode);
        for (int a = 0; a < k; ++a) {
            for (int b = 0; b < a; ++b) {
                final double mean = 0.5 * (hessian[a][b] + hessian[b][a]);
                hessian[a][b] = mean;
                hessian[b][a] = mean;
            }
        }
        return hessian;
    }

    /**
     * @return the lower Cholesky factor of a symmetric positive definite matrix
     * @throws IllegalStateException when the matrix is not positive definite
     */
    public static double[][] cholesky(double[][] matrix) {
        final int n = matrix.length;
        final double[][] lower = new double[n][n];
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j <= i; ++j) {
                double sum = matrix[i][j];
                for (int m = 0; m < j; ++m) {
                    sum -= lower[i][m] * lower[j][m];
                }
                if (i == j) {
                    if (!(sum > 0.0)) {
                        throw new IllegalStateException("The negative Hessian is not positive definite (pivot " + sum +
                                " at coordinate " + i + "): the point is not a maximum, or the density is flat there");
                    }
                    lower[i][i] = Math.sqrt(sum);
                } else {
                    lower[i][j] = sum / lower[j][j];
                }
            }
        }
        return lower;
    }

    // (L L^T)^-1 from the Cholesky factor L
    static double[][] inverse(double[][] lower) {
        final int n = lower.length;
        final double[][] inverseLower = new double[n][n]; // L^-1, lower triangular
        for (int j = 0; j < n; ++j) {
            inverseLower[j][j] = 1.0 / lower[j][j];
            for (int i = j + 1; i < n; ++i) {
                double sum = 0.0;
                for (int m = j; m < i; ++m) {
                    sum -= lower[i][m] * inverseLower[m][j];
                }
                inverseLower[i][j] = sum / lower[i][i];
            }
        }
        final double[][] inverse = new double[n][n];
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j <= i; ++j) {
                double sum = 0.0;
                for (int m = i; m < n; ++m) {
                    sum += inverseLower[m][i] * inverseLower[m][j];
                }
                inverse[i][j] = sum;
                inverse[j][i] = sum;
            }
        }
        return inverse;
    }

    private static double[][] copy(double[][] matrix) {
        final double[][] copy = new double[matrix.length][];
        for (int i = 0; i < matrix.length; ++i) {
            copy[i] = matrix[i].clone();
        }
        return copy;
    }
}

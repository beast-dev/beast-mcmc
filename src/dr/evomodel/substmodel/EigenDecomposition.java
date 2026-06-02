/*
 * EigenDecomposition.java
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

package dr.evomodel.substmodel;

import java.io.Serializable;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @Author Marc A. Suchard
 */
public class EigenDecomposition implements Serializable {

    public EigenDecomposition(double[] evec, double[] ievc, double[] eval) {
        Evec = evec;
        Ievc = ievc;
        Eval = eval;
    }

    public EigenDecomposition copy() {
        double[] evec = Evec.clone();
        double[] ievc = Ievc.clone();
        double[] eval = Eval.clone();

        return new EigenDecomposition(evec, ievc, eval);
    }

    public EigenDecomposition transpose() {
        // note: exchange e/ivec
        int dim = (int) Math.sqrt(Ievc.length);
        double[] evec = Ievc.clone();
        transposeInPlace(evec, dim);
        double[] ievc = Evec.clone();
        transposeInPlace(ievc, dim);
        double[] eval = Eval.clone();
        return new EigenDecomposition(evec, ievc, eval);
    }

    private static void transposeInPlace(double[] matrix, int n) {
        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                int index1 = i * n + j;
                int index2 = j * n + i;

                double temp = matrix[index1];
                matrix[index1] = matrix[index2];
                matrix[index2] = temp;
            }
        }
    }

    /**
     * This function returns the Eigen vectors.
     * @return the array
     */
    public final double[] getEigenVectors() {
        return Evec;
    }

    /**
     * This function returns the inverse Eigen vectors.
     * @return the array
     */
    public final double[] getInverseEigenVectors() {
        return Ievc;
    }

    /**
     * This function returns the Eigen values.
     * @return the Eigen values
     */
    public final double[] getEigenValues() {
        return Eval;
    }

    /**
     * This function returns the normalization factor
     * @return normalization factor
     */
    public final double getNormalization() { return normalization; }

    /**
     * This function rescales the eigen values; this is more stable than
     * rescaling the original Q matrix, also O(stateCount) instead of O(stateCount^2)
     */
    public void normalizeEigenValues(double scale) {
        this.normalization = scale;
        int dim = Eval.length;
        for (int i = 0; i < dim; i++) {
            Eval[i] /= scale;
        }
    }

    // Eigenvalues, eigenvectors, and inverse eigenvectors
    private final double[] Evec;
    private final double[] Ievc;
    private final double[] Eval;
    private double normalization = 1.0;
}

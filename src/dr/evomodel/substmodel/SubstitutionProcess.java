/*
 * SubstitutionProcess.java
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

package dr.evomodel.substmodel;

import dr.evolution.datatype.DataType;

/**
 * <b>model of sequence substitution (rate matrix)</b>.
 * provides a convenient interface for the computation of transition probabilities
 */
public interface SubstitutionProcess {

    /**
     * Get the complete transition probability matrix for the given distance.
     *
     * @param distance the time (branch length)
     * @param matrix   an array to store the matrix
     */
    void getTransitionProbabilities(double distance, double[] matrix);

    /**
     * This function returns the Eigen vectors.
     *
     * @return the array
     */
    EigenDecomposition getEigenDecomposition();

    /**
     * get the state frequencies
     *
     * @return the frequencies
     */
    FrequencyModel getFrequencyModel();

    /**
     * Get the infinitesimal rate matrix
     *
     * @return the rates
     */
    void getInfinitesimalMatrix(double[] matrix);

    /**
     * @return the data type
     */
    DataType getDataType();

    /**
     * @return if substitution model can return complex diagonalizations
     */
    boolean canReturnComplexDiagonalization();
}
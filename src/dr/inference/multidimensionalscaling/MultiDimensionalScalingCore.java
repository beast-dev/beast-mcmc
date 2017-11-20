/*
 * MultiDimensionalScalingCore.java
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

package dr.inference.multidimensionalscaling;

/**
 * MultiDimensionalScalingCore - An interface describing the core likelihood functions.
 *
 * @author Andrew Rambaut
 * @author Marc Suchard
 * @version $Id$
 */

public interface MultiDimensionalScalingCore {

    @SuppressWarnings("all")
    long USE_NATIVE_MDS = 1 << 0;       // 1
                                        // 2
    @SuppressWarnings("unused")
    long SINGLE_PRECISION = 1 << 2;     // 4
    @SuppressWarnings("unused")
    long MULTI_CORE = 1 << 3;           // 8
    @SuppressWarnings("unused")
    long OPENCL_VECTORIZATION = 1 << 4; // 16
    long LEFT_TRUNCATION = 1 << 5;      // 32

    /**
     * initializes arrays.
     */
    void initialize(int embeddingDimension, int locationCount, long flags);

    /**
     * sets the observation data
     */
    void setPairwiseData(double[] observations);

    /**
     *
     */
    void setParameters(double[] parameters);

    /**
     * Updates the location of a single entity
     */
    void updateLocation(int locationIndex, double[] location);

    /**
     * Calculates the log likelihood of the data given the current locations.
     */
    double calculateLogLikelihood();

    /**
     * Store current state
     */
    void storeState();

    /**
     * Restore the stored state
     */
    void restoreState();

    /**
     * Set all recalculate flags
     */
    void makeDirty();

    /**
     * Accept the proposed state
     */
    void acceptState();

    /**
     * Get gradient of MDS likelihood w.r.t. locations
     */
    void getGradient(double[] location);

    /**
     * Get pair-wise data
     */
    @SuppressWarnings("unused")
    double[] getPairwiseData();
}

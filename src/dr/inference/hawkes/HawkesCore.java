/*
 * HawkesCore.java
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

package dr.inference.hawkes;

/**
 * HawkesCore - An interface describing the core likelihood functions.
 *
 * @author Andrew Rambaut
 * @author Marc Suchard
 * @version $Id$
 */

public interface HawkesCore {

    @SuppressWarnings("all")
    long USE_NATIVE_HPH = 1 << 0;       // 1
                                        // 2
    @SuppressWarnings("unused")
    long SINGLE_PRECISION = 1 << 2;     // 4
    @SuppressWarnings("unused")
    long MULTI_CORE = 1 << 3;           // 8
    @SuppressWarnings("unused")
    long OPENCL_VECTORIZATION = 1 << 4; // 16

    /**
     * initializes arrays.
     */
    void initialize(int embeddingDimension, int locationCount, long flags);

    /**
     *
     */
    void setParameters(double[] parameters);

    /**
     *
     */
    void setTimesData(double[] timesData);

    /**
     *
     */
    void setRandomRates(double[] randomRates);

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
     * Get gradient of HPH likelihood w.r.t. locations
     */
    void getLocationGradient(double[] location);

    /**
     * Get gradient of HPH likelihood w.r.t. randomRates
     */
    void getRandomRatesGradient(double[] rate);

    int getInternalDimension();

    int DEFAULT_DEVICE = -1;

    class CoreInformation {
        public long flags = 0;
        public int deviceNumber = DEFAULT_DEVICE;
        public int numThreads = 1;
    }
}

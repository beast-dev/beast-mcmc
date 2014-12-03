/*
 * LikelihoodCore.java
 *
 * Copyright (c) 2002-2014 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.app.beagle.multidimensionalscaling;

/**
 * MultiDimensionalScalingCore - An interface describing the core likelihood functions.
 *
 * @author Andrew Rambaut
 * @version $Id$
 */

public interface MultiDimensionalScalingCore {

    /**
     * initializes arrays.
     */
    void initialize(int dimensionCount, int locationCount, int observationCount);

    /**
     * cleans up and deallocates arrays.
     */
    void finalize() throws Throwable;

    /**
     * sets the observation data
     */
    void setData(double[] observations, int[] observationTypes);

    /**
     * Allocates partials for a node
     */
    void updateLocations(int updateCount, double[] locations);

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
}

/*
 * MassivelyParallelMDSImpl.java
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
 * MassivelyParallelMDSImpl
 *
 * @author Marc Suchard
 * @author Andrew Rambaut
 * @version $Id$
 *          <p/>
 *          $HeadURL$
 *          <p/>
 *          $LastChangedBy$
 *          $LastChangedDate$
 *          $LastChangedRevision$
 */
public class MassivelyParallelMDSImpl implements MultiDimensionalScalingCore {

    @Override
    public void initialize(int embeddingDimension, int locationCount) {
        nativeInitialize(embeddingDimension, locationCount);
    }


    // embeddingDimension

    // dim rowColumnCount

    @Override
    public void setPairwiseData(double[] observations) {
        nativeSetPairwiseData(observations);
    }

    @Override
    public void setParameters(double[] parameters) {
        nativeSetParameters(parameters);
    }

    @Override
    public void updateLocation(int locationIndex, double[] location) {
        nativeUpdateLocation(locationIndex, location);
    }

    @Override
    public double calculateLogLikelihood() {
        return nativeCalculateLogLikelihood();
    }

    @Override
    public void storeState() {
        nativeStoreState();
    }

    @Override
    public void restoreState() {
        nativeRestoreState();
    }

    private native void nativeInitialize(int dimensionCount, int locationCount);

    private native void nativeSetPairwiseData(double[] observations);

    private native void nativeSetParameters(double[] parameters);

    private native void nativeUpdateLocation(int locationIndex, double[] location);

    private native double nativeCalculateLogLikelihood();

    private native void nativeStoreState();

    private native void nativeRestoreState();

}

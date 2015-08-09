/*
 * MassivelyParallelMDSImpl.java
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

    private NativeMDSSingleton singleton = null;
    private int instance = -1; // Get instance # via initialization
    private final long flags = 0;

    private static final long LEFT_TRUNCATION = 1 << 5;

    public MassivelyParallelMDSImpl() {
        singleton = NativeMDSSingleton.loadLibrary();
    }

    @Override
    public void initialize(int embeddingDimension, int locationCount, boolean isLeftTruncated) {
        long flags = this.flags;
        if (isLeftTruncated) {
            flags |= LEFT_TRUNCATION;
        }

        instance = singleton.initialize(embeddingDimension, locationCount, flags);
        this.observationCount = (locationCount * (locationCount - 1)) / 2;
        this.isLeftTruncated = isLeftTruncated;
    }

    @Override
    public void setPairwiseData(double[] observations) {
        singleton.setPairwiseData(instance, observations);
    }

    @Override
    public void setParameters(double[] parameters) {
        precision = parameters[0];
        singleton.setParameters(instance, parameters); // Necessary for truncation
    }

    @Override
    public void updateLocation(int locationIndex, double[] location) {
        singleton.updateLocations(instance, locationIndex, location);
    }

    @Override
    public double calculateLogLikelihood() {
        double sumOfSquaredResiduals = singleton.getSumOfSquaredResiduals(instance);

        double logLikelihood = (0.5 * Math.log(precision) * observationCount) -
                        (0.5 * precision * sumOfSquaredResiduals);

        if (isLeftTruncated) {
            logLikelihood -= singleton.getSumOfLogTruncations(instance);
        }

        return logLikelihood;
    }

    @Override
    public void storeState() {
        singleton.storeState(instance);
        storedPrecision = precision;
    }

    @Override
    public void restoreState() {
        singleton.restoreState(instance);
        precision = storedPrecision;
    }

    @Override
    public void makeDirty() {
        singleton.makeDirty(instance);
    }

    private int observationCount;
    private double precision;
    private double storedPrecision;
    private boolean isLeftTruncated;

}

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

package dr.inference.multidimensionalscaling;

import dr.math.matrixAlgebra.Vector;

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

    MassivelyParallelMDSImpl() {
        singleton = NativeMDSSingleton.loadLibrary();
    }

    @Override
    public void initialize(int embeddingDimension, int locationCount, long flags) {
        instance = singleton.initialize(embeddingDimension, locationCount, flags);
        this.observationCount = (locationCount * (locationCount - 1)) / 2;
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
    public double[] getPairwiseData() {
        return singleton.getPairwiseData(instance);
    }

    @Override
    public void updateLocation(int locationIndex, double[] location) {
        singleton.updateLocations(instance, locationIndex, location);
    }

    @Override
    public double calculateLogLikelihood() {
        double sumOfIncrements = singleton.getSumOfIncrements(instance);

        return 0.5 * (Math.log(precision) - Math.log(2 * Math.PI)) * observationCount - sumOfIncrements;
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
    public void acceptState() {
        singleton.acceptState(instance);
    }

    @Override
    public void getGradient(double[] location) {
        singleton.getLocationGradient(instance, location);

        if (CHECK_GRADIENT) {
            for (double x : location) {
                if (Double.isNaN(x) || Double.isInfinite(x)) {
                    System.err.println("Poor gradient value: " + x);
                    System.err.println(new Vector(location));
                    if (CHECK_GRADIENT_KILL) {
                        System.exit(-1);
                    } else {
                        break;
                    }
                }
            }
        }
    }

    @Override
    public void makeDirty() {
        singleton.makeDirty(instance);
    }

    private int observationCount;
    private double precision;
    private double storedPrecision;

    private static final boolean CHECK_GRADIENT = false;
    private static final boolean CHECK_GRADIENT_KILL = true;
}

/*
 * MassivelyParallelHawkesImpl.java
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

import dr.math.matrixAlgebra.Vector;

import static dr.inference.hawkes.NativeHPHSingleton.HPH_RESOURCE;
import static dr.inference.hawkes.NativeHPHSingleton.THREADS;

/**
 * MassivelyParallelHPHImpl
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
public class MassivelyParallelHPHImpl implements HawkesCore {

    private NativeHPHSingleton singleton;
    private int instance = -1; // Get instance # via initialization

    private final CoreInformation information;

    MassivelyParallelHPHImpl() {
        singleton = NativeHPHSingleton.loadLibrary();

        information = new CoreInformation();

        String resource = System.getProperty(HPH_RESOURCE);
        if (resource != null) {
            try {
                int number = Integer.parseInt(resource);
                if (number > 0) {
                    information.deviceNumber = number - 1;
                }
            } catch (NumberFormatException exception) {
                throw new RuntimeException("Unable to parse '" + HPH_RESOURCE + "' environmental property");
            }
        }

        String r = System.getProperty(THREADS);
        if (r != null) {
            try {
                information.numThreads = Integer.parseInt(r.trim());
            } catch (NumberFormatException exception) {
                throw new RuntimeException("Unable to parse '" + THREADS + "' environmental property");
            }
        }
    }

    @Override
    public void initialize(int embeddingDimension, int locationCount, long flags) {
        information.flags = flags;
        instance = singleton.initialize(embeddingDimension, locationCount, information);
        this.observationCount = (locationCount * (locationCount - 1)) / 2;
    }

    @Override
    public void setTimesData(double[] timesData) {
        singleton.setTimesData(instance, timesData);
    }

    @Override
    public void setRandomRates(double[] randomRates) {
        singleton.setRandomRates(instance, randomRates);
    }

    @Override
    public void setParameters(double[] parameters) {
        sigmaXprec = parameters[0];
        tauXprec = parameters[1];
        tauTprec = parameters[2];
        omega = parameters[3];
        theta = parameters[4];
        mu0 = parameters[5];

        singleton.setParameters(instance, parameters);
    }

    @Override
    public void updateLocation(int locationIndex, double[] location) {
        singleton.updateLocations(instance, locationIndex, location);
    }

    @Override
    public double calculateLogLikelihood() {
        double sumOfIncrements = singleton.getSumOfIncrements(instance);

        return sumOfIncrements;
    }

    @Override
    public void storeState() {
        singleton.storeState(instance);
        storedSigmaXprec = sigmaXprec;
        storedTauXprec = tauXprec;
        storedTauTprec = tauTprec;
        storedOmega = omega;
        storedTheta = theta;
        storedMu0 = mu0;    }

    @Override
    public void restoreState() {
//        singleton.restoreState(instance);
        sigmaXprec = storedSigmaXprec;
        tauXprec = storedTauXprec;
        tauTprec = storedTauTprec;
        omega = storedOmega;
        theta = storedTheta;
        mu0 = storedMu0;
    }

    @Override
    public void acceptState() {
        singleton.acceptState(instance);
    }

    @Override
    public void getLocationGradient(double[] location) {
        singleton.getLocationGradient(instance, location);

        if (CHECK_GRADIENT) {
            checkGradient(location);
        }
    }

    @Override
    public void getRandomRatesGradient(double[] rate) {
        singleton.getRandomRatesGradient(instance, rate);

        if (CHECK_GRADIENT) {
            checkGradient(rate);
        }
    }

    @Override
    public void getRandomRatesHessian(double[] rate) {
        singleton.getRandomRatesHessian(instance, rate);

        if (CHECK_GRADIENT) {
            checkGradient(rate);
        }
    }

    private void checkGradient(double[] array) {
        for (double x : array) {
                if (Double.isNaN(x) || Double.isInfinite(x)) {
                    System.err.println("Poor gradient value: " + x);
                    System.err.println(new Vector(array));
                    if (CHECK_GRADIENT_KILL) {
                        System.exit(-1);
                    } else {
                        break;
                    }
                }
            }
    }

    @Override
    public void makeDirty() {
//        singleton.makeDirty(instance);
    }

    @Override
    public int getInternalDimension() { return singleton.getInternalDimension(instance); }

    private int observationCount;
    private double tauXprec;
    private double storedTauXprec;
    private double tauTprec;
    private double storedTauTprec;
    private double sigmaXprec;
    private double storedSigmaXprec;
    private double omega;
    private double storedOmega;
    private double theta;
    private double storedTheta;
    private double mu0;
    private double storedMu0;


    private static final boolean CHECK_GRADIENT = false;
    private static final boolean CHECK_GRADIENT_KILL = true;
}

/*
 * MultiDimensionalScalingCoreImpl.java
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

import dr.math.distributions.NormalDistribution;

/**
 * MultiDimensionalScalingCoreImpl
 *
 * @author Andrew Rambaut
 * @author Marc Suchard
 * @version $Id$
 *
 * $HeadURL$
 *
 * $LastChangedBy$
 * $LastChangedDate$
 * $LastChangedRevision$
 */

public class MultiDimensionalScalingCoreImpl implements MultiDimensionalScalingCore {
    @Override
    public void initialize(int embeddingDimension, int locationCount) {
        this.embeddingDimension = embeddingDimension;
        this.locationCount = locationCount;
        this.observationCount = locationCount * locationCount;

        observations = new double[locationCount * locationCount];
        locationUpdated = new boolean[locationCount];

        locations = new double[locationCount][embeddingDimension];
    }

    @Override
    public void setPairwiseData(double[] observations) {
        if (observations.length != observationCount) {
            throw new RuntimeException("Observation data is not the correct dimension");
        }
        System.arraycopy(observations, 0, this.observations, 0, observationCount);
    }

    @Override
    public void setParameters(double[] parameters) {
        precision = parameters[0];
    }

    @Override
    public void updateLocation(int locationIndex, double[] location) {
        locationUpdated[locationIndex] = true;
        likelihoodKnown = false;
        distancesKnown = false;
        residualsKnown = false;

        if (location.length != embeddingDimension) {
            throw new RuntimeException("Location is not the correct dimension");
        }
        System.arraycopy(location, 0, locations[locationIndex], 0, embeddingDimension);
    }

    @Override
    public double calculateLogLikelihood() {
        if (!likelihoodKnown) {
            if (locationsUpdated) {
                calculateDistances();
                residualsKnown = false;
            }

            double precision = mdsPrecisionParameter.getParameterValue(0);

            if (!residualsKnown) {
                sumOfSquaredResiduals = calculateSumOfSquaredResiduals();
            }

            // totalNonMissingCount should be totalObservedCount (not > or < threshold)
            double logLikelihood = (0.5 * Math.log(precision) * pointObservationCount) - (0.5 * precision * sumOfSquaredResiduals);

            if (thresholdCount > 0) {
                if (!thresholdsKnown) {
                    thresholdSum = calculateThresholdObservations(precision);
                }
                logLikelihood += thresholdSum;
            }

            if (isLeftTruncated) {
                if (!truncationsKnown) {
                    calculateTruncations(precision);
                }
                truncationSum = calculateTruncationSum();
                logLikelihood -= truncationSum;
            }

            likelihoodKnown = true;

            for (int i = 0; i < locationUpdated.length; i++) {
                locationUpdated[i] = false;
            }
        }

        return logLikelihood;
    }

    @Override
    public void storeState() {

    }

    @Override
    public void restoreState() {

    }

    protected void calculateTruncations(double precision) {
        double sd = 1.0 / Math.sqrt(precision);
        for (int i = 0; i < distanceCount; i++) {
            if (distanceUpdated[i]) {
                truncations[i] = NormalDistribution.cdf(distances[i], 0.0, sd, true);
            }
        }
        truncationsKnown = true;
    }

    protected double calculateTruncationSum() {
        double sum = 0.0;
        for (int i = 0; i < observationCount; i++) {
            int dist = getDistanceIndexForObservation(i);
            if (dist != -1) {
                sum += truncations[dist];
            } else {
                sum += Math.log(0.5);
            }
        }
        return sum;
    }

    protected double calculateSumOfSquaredResiduals() {
        double sum = 0.0;
        for (int i = 0; i < observationCount; i++) {
            if (observationTypes[i] == ObservationType.POINT) {
                // Only increment sum if dataTable[i][j] is observed (not > or < threshold)
                double residual;
                int dist = getDistanceIndexForObservation(i);
                if (dist == -1) {
                    // -1 denotes a distance to self (i.e., 0)
                    residual = - observations[i];
                } else {
                    residual = distances[dist] - observations[i];
                }
                sum += residual * residual;
            }
        }
        residualsKnown = true;
        return sum;
    }

    protected void calculateDistances() {
        int k = 0;
        for (int x = 0; x < locationCount; x++) {
            for (int y = x + 1; y < locationCount; y++) {
                // the diagonal (x=y) is always zero so don't update it
                if (locationUpdated[x] || locationUpdated[y]) {
                    distances[k] = calculateDistance(
                            locationsParameter.getParameter(x),
                            locationsParameter.getParameter(y));
                    distanceUpdated[k] = true;
                }
                k++;
            }
        }
        distancesKnown = true;
    }

    private int embeddingDimension;
    private boolean isLeftTruncated = false;
    private int locationCount;
    private int observationCount;
    private double precision;

    private double[] observations;
    private double[][] locations;

    private boolean[] locationUpdated;

    private boolean likelihoodKnown = false;
    private double logLikelihood;
    private double storedLogLikelihood;

    protected boolean distancesKnown = false;
    private double sumOfSquaredResiduals;
    private double storedSumOfSquaredResiduals;

    protected boolean residualsKnown = false;

    protected boolean truncationsKnown = false;
    private double truncationSum;
    private double storedTruncationSum;
    private double[] truncations;
    private double[] storedTruncations;

}

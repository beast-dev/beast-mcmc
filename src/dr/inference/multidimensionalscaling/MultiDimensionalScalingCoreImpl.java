/*
 * MultiDimensionalScalingCoreImpl2.java
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
    public void initialize(int embeddingDimension, int locationCount, long flags) {
        this.embeddingDimension = embeddingDimension;
        this.locationCount = locationCount;
        this.observationCount = (locationCount * (locationCount - 1)) / 2;

        observations = new double[locationCount][locationCount];
        increments = new double[locationCount][locationCount];
        storedIncrements = null;
        incrementsKnown = false;
        sumOfIncrementsKnown = false;

        isLeftTruncated = (flags & MultiDimensionalScalingCore.LEFT_TRUNCATION) != 0;

        updatedLocation = -1;

        locations = new double[locationCount][embeddingDimension];
        storedLocations = new double[locationCount][embeddingDimension];
    }

    @Override
    public void setPairwiseData(double[] observations) {
        if (observations.length != (locationCount * locationCount)) {
            throw new RuntimeException("Observation data is not the correct dimension");
        }

        int k = 0;
        for (int i = 0; i < locationCount; i++) {
            System.arraycopy(observations, k, this.observations[i], 0, locationCount);
            k += locationCount;
        }
    }

    @Override
    public double[] getPairwiseData() {
        double[] data = new double[locationCount * locationCount];
        int k = 0;
        for (int i = 0; i < locationCount; ++i) {
            System.arraycopy(observations[i], 0, data, k, locationCount);
            k += locationCount;
        }
        return data;
    }

    @Override
    public void setParameters(double[] parameters) {
        precision = parameters[0];

        // Handle truncations
        if (isLeftTruncated) {
            incrementsKnown = false;
            sumOfIncrementsKnown = false;
        }
    }

    @Override
    public void updateLocation(int locationIndex, double[] location) {
        if (updatedLocation != -1 || locationIndex == -1) {
            // more than one location updated - do a full re-computation
            incrementsKnown = false;
            storedIncrements = null;
         }

        if (locationIndex != -1) {
            updatedLocation = locationIndex;

            if (location.length != embeddingDimension) {
                throw new RuntimeException("Location is not the correct dimension");
            }

            System.arraycopy(location, 0, locations[locationIndex], 0, embeddingDimension);

        } else {
            if (location.length != embeddingDimension * locationCount) {
                throw new RuntimeException("Location is the not correct dimension");
            }

            int offset = 0;
            for (int i = 0; i < locationCount; ++i) {
                System.arraycopy(location, offset, locations[i], 0, embeddingDimension);
                offset += embeddingDimension;
            }
        }

        sumOfIncrementsKnown = false;
    }

    @Override
    public double calculateLogLikelihood() {
        if (!sumOfIncrementsKnown) {

            if (!incrementsKnown) {
                computeSumOfSquaredResiduals();
            } else {
                updateSumOfSquaredResiduals();
            }
            sumOfIncrementsKnown = true;
        }

        double logLikelihood = 0.5 * (Math.log(precision) - Math.log(2 * Math.PI)) * observationCount;

        if (isLeftTruncated) {
            logLikelihood -= sumOfIncrements; // If truncated, then values on difference scale
        } else {
            logLikelihood -= 0.5 * precision * sumOfIncrements;
        }

        return logLikelihood;
    }

    @Override
    public void storeState() {
        // Handle residuals
        storedSumOfIncrements = sumOfIncrements;
        storedIncrements = null;

        // Handle locations
        for (int i = 0; i < locationCount; i++) {
            System.arraycopy(locations[i], 0 , storedLocations[i], 0, embeddingDimension);
        }
        updatedLocation = -1;

        // Handle precision
        storedPrecision = precision;
    }

    @Override
    public void restoreState() {
        // Handle residuals
        sumOfIncrements = storedSumOfIncrements;
        sumOfIncrementsKnown = true;

        if (storedIncrements != null) {
            System.arraycopy(storedIncrements, 0 , increments[updatedLocation], 0, locationCount);
//            for (int j = 0; j < locationCount; j++) { // Do not write transposed values
//                increments[j][updatedLocation] = storedIncrements[j];
//            }
            incrementsKnown = true;
        } else {
            incrementsKnown = false;
        }

        // Handle locations
        double[][] tmp1 = storedLocations;
        storedLocations = locations;
        locations = tmp1;

        // Handle precision
        precision = storedPrecision;
    }

    @Override
    public void acceptState() {
        if (storedIncrements != null) {
            for (int j = 0; j < locationCount; ++j) {
                increments[j][updatedLocation] = increments[updatedLocation][j];
            }
        }
    }

    @Override
    public void getGradient(double[] location) {
        throw new RuntimeException("Not yet implemented.");
    }

    @Override
    public void makeDirty() {
        sumOfIncrementsKnown = false;
        incrementsKnown = false;
    }

    private void computeSumOfSquaredResiduals() {

        final double oneOverSd = Math.sqrt(precision);
        final double scale = 0.5 * precision;

        // OLD
        sumOfIncrements = 0.0;
        for (int i = 0; i < locationCount; i++) {

            for (int j = 0; j < locationCount; j++) {
                double distance = calculateDistance(locations[i], locations[j]);
                double residual = distance - observations[i][j];
                double increment = residual * residual;
                if (isLeftTruncated) {
                    increment = scale * increment;
                    if (i != j) {
                        increment += computeTruncation(distance, oneOverSd);
//                        increment += computeTruncation(Math.sqrt(residual * residual), precision, oneOverSd); // OLD .. believed incorrect
                    }
                }
                increments[i][j] = increment;
//                increments[j][i] = increment; // Do not write transposed values
                sumOfIncrements += increment;
            }
        }

        sumOfIncrements /= 2;

        incrementsKnown = true;
        sumOfIncrementsKnown = true;
    }

    private void updateSumOfSquaredResiduals() {

        final double oneOverSd = Math.sqrt(precision);
        final double scale = 0.5 * precision;

        double delta = 0.0;

        final int i = updatedLocation;

        storedIncrements = new double[locationCount];
        System.arraycopy(increments[i], 0, storedIncrements, 0, locationCount);

        for (int j = 0; j < locationCount; j++) {
            double distance = calculateDistance(locations[i], locations[j]);
            double residual = distance - observations[i][j];
            double increment = residual * residual;

            if (isLeftTruncated) {
                increment = scale * increment;
                if (i != j) {
                    increment += computeTruncation(distance, oneOverSd);
//                    increment += computeTruncation(Math.sqrt(residual * residual), precision, oneOverSd); // OLD .. believed incorrect
                }
            }

            delta += increment - increments[i][j];
            increments[i][j] = increment;
//            increments[j][i] = increment; // Do not write transposed values
        }

        sumOfIncrements += delta;
    }

    private double calculateDistance(double[] X, double[] Y) {
        double sum = 0.0;
        for (int i = 0; i < embeddingDimension; i++) {
            double difference = X[i] - Y[i];
            sum += difference * difference;
        }
        return Math.sqrt(sum);
    }

    private double computeTruncation(double mean, double oneOverSd) {
        return NormalDistribution.standardCDF(mean * oneOverSd, true); // Should be standardCDF(mean / sd, true);
    }

    private int embeddingDimension;
    private boolean isLeftTruncated = false;
    private int locationCount;
    private int observationCount;
    private double precision;
    private double storedPrecision;

    private int updatedLocation = -1;

    private double[][] observations;
    private double[][] locations;
    private double[][] storedLocations;

    private boolean incrementsKnown = false;

    private boolean sumOfIncrementsKnown = false;
    private double[][] increments;

    private double[] storedIncrements;

    private double sumOfIncrements;
    private double storedSumOfIncrements;

}

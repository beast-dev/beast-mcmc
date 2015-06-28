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

public class MultiDimensionalScalingCoreImpl2 implements MultiDimensionalScalingCore {

    @Override
    public void initialize(int embeddingDimension, int locationCount, boolean isLeftTruncated) {
        this.embeddingDimension = embeddingDimension;
        this.locationCount = locationCount;
        this.observationCount = (locationCount * (locationCount - 1)) / 2;
        this.isLeftTruncated = isLeftTruncated;

        observations = new double[locationCount][locationCount];
        squaredResiduals = new double[locationCount][locationCount];
        storedSquaredResiduals = null;
        residualsKnown = false;
        sumOfSquaredResidualsKnown = false;

        if (isLeftTruncated) {
            truncations = new double[locationCount][locationCount];
            storedTruncations = null;
            truncationsKnown = false;
            sumOfTruncationsKnown = false;
        }

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
    public void setParameters(double[] parameters) {
        precision = parameters[0];

        // Handle truncations
        truncationsKnown = false;
        sumOfTruncationsKnown = false;
    }

    @Override
    public void updateLocation(int locationIndex, double[] location) {
        if (updatedLocation != -1 || locationIndex == -1) {
            // more than one location updated - do a full recomputation
            residualsKnown = false;
            storedSquaredResiduals = null;

            // Handle truncations
            truncationsKnown = false;
            storedTruncations = null;
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

        sumOfSquaredResidualsKnown = false;

        // Handle truncation
        sumOfTruncationsKnown = false;
    }

    @Override
    public double calculateLogLikelihood() {
        if (!sumOfSquaredResidualsKnown) {

            if (!residualsKnown) {
                computeSumOfSquaredResiduals();
            } else {
                updateSumOfSquaredResiduals();
                if (REPORT_ROUNDOFF) {
                    // Report round-off error
                    double storedSumOfSquaredResults = sumOfSquaredResiduals;
                    computeSumOfSquaredResiduals();
                    if (Math.abs(storedSumOfSquaredResults - sumOfSquaredResiduals) > 1E-6) {
                        System.err.println(storedSumOfSquaredResults);
                        System.err.println(sumOfSquaredResiduals);
                        System.err.println(storedSumOfSquaredResults - sumOfSquaredResiduals);
                        System.err.println("");
                    }
                }
            }
            sumOfSquaredResidualsKnown = true;
        }

        double logLikelihood = (0.5 * Math.log(precision) * observationCount) -
                (0.5 * precision * sumOfSquaredResiduals);

        if (isLeftTruncated) {
            if (!sumOfTruncationsKnown) {

                if (!truncationsKnown) {
                    computeSumOfTruncations();
                } else {
                    updateSumOfTruncations();
                }
                sumOfTruncationsKnown = true;
            }
            logLikelihood -= truncationSum;
        }

        return logLikelihood;
    }

    @Override
    public void storeState() {
        // Handle residuals
        storedSumOfSquaredResiduals = sumOfSquaredResiduals;
        storedSquaredResiduals = null;

        // Handle locations
        for (int i = 0; i < locationCount; i++) {
            System.arraycopy(locations[i], 0 , storedLocations[i], 0, embeddingDimension);
        }
        updatedLocation = -1;

        // Handle precision
        storedPrecision = precision;

        // Handle truncations
        if (isLeftTruncated) {
            storedTruncationSum = truncationSum;
            storedTruncations = null;
        }
    }

    @Override
    public void restoreState() {
        // Handle residuals
        sumOfSquaredResiduals = storedSumOfSquaredResiduals;
        sumOfSquaredResidualsKnown = true;

        if (storedSquaredResiduals != null) {
            System.arraycopy(storedSquaredResiduals, 0 , squaredResiduals[updatedLocation], 0, locationCount);
//            for (int j = 0; j < locationCount; j++) { // Do not write transposed values
//                squaredResiduals[j][updatedLocation] = storedSquaredResiduals[j];
//            }
            residualsKnown = true;
        } else {
            residualsKnown = false;
        }

        // Handle locations
        double[][] tmp1 = storedLocations;
        storedLocations = locations;
        locations = tmp1;

        // Handle precision
        precision = storedPrecision;

        // Handle truncations
        if (isLeftTruncated) {
            truncationSum = storedTruncationSum;
            sumOfTruncationsKnown = true;

            if (storedTruncations != null) {
                System.arraycopy(storedTruncations, 0, truncations[updatedLocation], 0, locationCount);
//                for (int j = 0; j < locationCount; ++j) { // Do not write transposed values
//                    truncations[j][updatedLocation] = storedTruncations[j];
//                }
                truncationsKnown = true;
            } else {
                truncationsKnown = false;
            }
        }
    }

    @Override
    public void acceptState() {
        if (storedSquaredResiduals != null) {
            for (int j = 0; j < locationCount; ++j) {
                squaredResiduals[j][updatedLocation] = squaredResiduals[updatedLocation][j];
            }
        }

        if (isLeftTruncated) {
            if (storedTruncations != null) {
                for (int j = 0; j < locationCount; ++j) { // Do not write transposed values
                    truncations[j][updatedLocation] = truncations[updatedLocation][j];
                }
            }
        }
    }

    @Override
    public void makeDirty() {
        sumOfSquaredResidualsKnown = false;
        residualsKnown = false;

        sumOfTruncationsKnown = false;
        truncationsKnown = false;
    }

    protected void computeSumOfSquaredResiduals() {

        // OLD
        sumOfSquaredResiduals = 0.0;
        for (int i = 0; i < locationCount; i++) {

            for (int j = 0; j < locationCount; j++) {
                double distance = calculateDistance(locations[i], locations[j]);
                double residual = distance - observations[i][j];
                double squaredResidual = residual * residual;
                squaredResiduals[i][j] = squaredResidual;
//                squaredResiduals[j][i] = squaredResidual; // Do not write transposed values
                sumOfSquaredResiduals += squaredResidual;
            }
        }

        sumOfSquaredResiduals /= 2;

        // New   TODO
//        sumOfSquaredResiduals = 0.0;
//         for (int i = 0; i < locationCount; i++) {
//
//             for (int j = i + 1; j < locationCount; j++) {
//                 double distance = calculateDistance(locations[i], locations[j]);
//                 double residual = distance - observations[i][j];
//                 double squaredResidual = residual * residual;
//                 squaredResiduals[i][j] = squaredResidual;
//                 squaredResiduals[j][i] = squaredResidual;
//                 sumOfSquaredResiduals += squaredResidual;
//             }
//         }

        residualsKnown = true;
        sumOfSquaredResidualsKnown = true;
    }

    protected void computeSumOfTruncations() {

        final double oneOverSd = Math.sqrt(precision);

        truncationSum = 0.0;
        for (int i = 0; i < locationCount; i++) {

            for (int j = 0; j < locationCount; j++) {
                double squaredResidual = squaredResiduals[i][j]; // Note just written above, save transaction
                double truncation = (i == j) ? 0.0 : computeTruncation(squaredResidual, precision, oneOverSd);
                truncations[i][j] =  truncation;
//                truncations[j][i] = truncation;
                truncationSum += truncation;
            }
        }

        truncationSum /= 2;

        truncationsKnown = true;
        sumOfTruncationsKnown = true;
    }

    protected void updateSumOfSquaredResiduals() {
        double delta = 0.0;

        int i = updatedLocation;

        storedSquaredResiduals = new double[locationCount];
        System.arraycopy(squaredResiduals[i], 0, storedSquaredResiduals, 0, locationCount);

        for (int j = 0; j < locationCount; j++) {
            double distance = calculateDistance(locations[i], locations[j]);
            double residual = distance - observations[i][j];
            double squaredResidual = residual * residual;

            delta += squaredResidual - squaredResiduals[i][j];

            squaredResiduals[i][j] = squaredResidual;
//            squaredResiduals[j][i] = squaredResidual; // Do not write transposed values
        }

        sumOfSquaredResiduals += delta;
    }

    protected void updateSumOfTruncations() {
        final double oneOverSd = Math.sqrt(precision);
        double delta = 0.0;

        int i = updatedLocation;

        storedTruncations = new double[locationCount];
        System.arraycopy(truncations[i], 0, storedTruncations, 0, locationCount);

        for (int j = 0; j < locationCount; j++) {

            double squaredResidual = squaredResiduals[i][j];
            double truncation = (i == j) ? 0.0 : computeTruncation(squaredResidual, precision, oneOverSd);

            delta += truncation - truncations[i][j];

            truncations[i][j] = truncation;
//            truncations[j][i] = truncation; // Do not write transposed values
        }

        truncationSum += delta;
    }

    protected double calculateDistance(double[] X, double[] Y) {
        double sum = 0.0;
        for (int i = 0; i < embeddingDimension; i++) {
            double difference = X[i] - Y[i];
            sum += difference * difference;
        }
        return Math.sqrt(sum);
    }

    protected double computeTruncation(double squaredResidual, double precision, double oneOverSd) {
        return NormalDistribution.standardCDF(Math.sqrt(squaredResidual) * oneOverSd, true);
    }

//    protected void calculateTruncations(double precision) {
//        double sd = 1.0 / Math.sqrt(precision);
//        for (int i = 0; i < distanceCount; i++) {
//            if (distanceUpdated[i]) {
//                truncations[i] = NormalDistribution.cdf(distances[i], 0.0, sd, true);
//            }
//        }
//        truncationsKnown = true;
//    }
//
//    protected double calculateTruncationSum() {
//        double sum = 0.0;
//        for (int i = 0; i < observationCount; i++) {
//            int dist = getDistanceIndexForObservation(i);
//            if (dist != -1) {
//                sum += truncations[dist];
//            } else {
//                sum += Math.log(0.5);
//            }
//        }
//        return sum;
//    }

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

    private boolean residualsKnown = false;

    private boolean sumOfSquaredResidualsKnown = false;
    private double[][] squaredResiduals;

    private double[] storedSquaredResiduals;

    private double sumOfSquaredResiduals;
    private double storedSumOfSquaredResiduals;

    private boolean truncationsKnown = false;
    private boolean sumOfTruncationsKnown = false;

    private double truncationSum;
    private double storedTruncationSum;
    private double[][] truncations;
    private double[] storedTruncations;

    private static boolean REPORT_ROUNDOFF = false;

}

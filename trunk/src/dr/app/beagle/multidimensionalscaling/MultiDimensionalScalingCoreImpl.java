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
        this.observationCount = (locationCount * (locationCount - 1)) / 2;

        observations = new double[locationCount][locationCount];
        squaredResiduals = new double[locationCount][locationCount];
        storedSquaredResiduals = new double[locationCount][locationCount];
        residualsKnown = false;
        sumOfSquaredResidualsKnown = false;

        locationUpdated = new boolean[locationCount];
        for (int i = 0; i < locationUpdated.length; i++) {
            locationUpdated[i] = true;
        }

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
    }

    @Override
    public void updateLocation(int locationIndex, double[] location) {
        if (location.length != embeddingDimension) {
            throw new RuntimeException("Location is not the correct dimension");
        }

        System.arraycopy(location, 0, locations[locationIndex], 0, embeddingDimension);
        locationUpdated[locationIndex] = true;
        sumOfSquaredResidualsKnown = false;
    }

    @Override
    public double calculateLogLikelihood() {
        if (!sumOfSquaredResidualsKnown) {
            if (!residualsKnown) {
                computeSumOfSquaredResiduals();
            } else {
                updateSumOfSquaredResiduals();
            }
            sumOfSquaredResidualsKnown = true;
        }

        double logLikelihood = (0.5 * Math.log(precision) * observationCount) -
                (0.5 * precision * sumOfSquaredResiduals);

        if (isLeftTruncated) {
            throw new UnsupportedOperationException("Truncations not implemented");
//                if (!truncationsKnown) {
//                    calculateTruncations(precision);
//                }
//                truncationSum = calculateTruncationSum();
//                logLikelihood -= truncationSum;
        }

        for (int i = 0; i < locationUpdated.length; i++) {
            locationUpdated[i] = false;
        }

        return logLikelihood;
    }

    @Override
    public void storeState() {
        storedSumOfSquaredResiduals = sumOfSquaredResiduals;
        for (int i = 0; i < locationCount; i++) {
            System.arraycopy(squaredResiduals[i], 0 , storedSquaredResiduals[i], 0, locationCount);
        }

        for (int i = 0; i < locationCount; i++) {
            System.arraycopy(locations[i], 0 , storedLocations[i], 0, embeddingDimension);
        }

        storedPrecision = precision;
    }

    @Override
    public void restoreState() {
        sumOfSquaredResiduals = storedSumOfSquaredResiduals;
        sumOfSquaredResidualsKnown = true;

        double[][] tmp = storedSquaredResiduals;
        storedSquaredResiduals = squaredResiduals;
        squaredResiduals = tmp;

        tmp = storedLocations;
        storedLocations = locations;
        locations = tmp;

        precision = storedPrecision;

        residualsKnown = true;
    }

    @Override
    public void makeDirty() {
        sumOfSquaredResidualsKnown = false;
        residualsKnown = false;


    }

    protected void computeSumOfSquaredResiduals() {
        sumOfSquaredResiduals = 0.0;
        for (int i = 0; i < locationCount; i++) {
            for (int j = i + 1; j < locationCount; j++) {
                double distance = calculateDistance(locations[i], locations[j]);
                double residual = distance - observations[i][j];
                double squaredResidual = residual * residual;
                squaredResiduals[i][j] = squaredResidual;
                squaredResiduals[j][i] = squaredResidual;
                sumOfSquaredResiduals += squaredResidual;
            }
        }

        residualsKnown = true;
        sumOfSquaredResidualsKnown = true;
    }

    protected void updateSumOfSquaredResiduals() {
        double delta = 0.0;

        double[] oldSquaredResidualRow = new double[locationCount];

        for (int i = 0; i < locationCount; i++) {
            if (locationUpdated[i]) {

                System.arraycopy(squaredResiduals[i], 0, oldSquaredResidualRow, 0, locationCount);

                // if location i is updated, calculate the residuals to all js
                // also sum the change in sum residual
                for (int j = 0; j < locationCount; j++) {
                    double distance = calculateDistance(locations[i], locations[j]);
                    double residual = distance - observations[i][j];
                    double squaredResidual = residual * residual;

                    delta += squaredResidual - oldSquaredResidualRow[j];

                    squaredResiduals[i][j] = squaredResidual;
                    squaredResiduals[j][i] = squaredResidual;
                }
            }
        }

        sumOfSquaredResiduals += delta;
    }

    protected double calculateDistance(double[] X, double[] Y) {
        double sum = 0.0;
        for (int i = 0; i < embeddingDimension; i++) {
            double difference = X[i] - Y[i];
            sum += difference * difference;
        }
        return Math.sqrt(sum);
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


    private double[][] observations;
    private double[][] locations;
    private double[][] storedLocations;

    private boolean[] locationUpdated;

    private boolean residualsKnown = false;

    private boolean sumOfSquaredResidualsKnown = false;
    private double[][] squaredResiduals;
    private double[][] storedSquaredResiduals;
    private double sumOfSquaredResiduals;
    private double storedSumOfSquaredResiduals;

    private boolean truncationsKnown = false;
    private double truncationSum;
    private double storedTruncationSum;
    private double[] truncations;
    private double[] storedTruncations;

}

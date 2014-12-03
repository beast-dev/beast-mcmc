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
        this.observationCount = (locationCount * (locationCount - 1)) / 2;

        observations = new double[locationCount][locationCount];
        squaredResiduals = new double[2][locationCount][locationCount];
        residualFlag = 0;
        residualsKnown = false;
        sumOfSquaredResidualsKnown = false;

        locationUpdated = new boolean[locationCount];
        for (int i = 0; i < locationUpdated.length; i++) {
            locationUpdated[i] = true;
        }

        locations = new double[locationCount][embeddingDimension];
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
        makeDirty();

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
        storedResidualFlag = residualFlag;
        storedSumOfSquaredResiduals = sumOfSquaredResiduals;
    }

    @Override
    public void restoreState() {
        residualFlag = storedResidualFlag;
        sumOfSquaredResiduals = storedSumOfSquaredResiduals;

        residualsKnown = true;
    }

    public void makeDirty() {
        sumOfSquaredResidualsKnown = false;
        residualsKnown = false;


    }

    protected void computeSumOfSquaredResiduals() {
        for (int i = 0; i < locationCount; i++) {
            for (int j = i + 1; j < locationCount; j++) {
                double distance = calculateDistance(locations[i], locations[j]);
                double residual = distance - observations[i][j];
                double squaredResidual = residual * residual;
                squaredResiduals[0][i][j] = squaredResidual;
                squaredResiduals[1][i][j] = squaredResidual;
                sumOfSquaredResiduals += squaredResidual;
            }
        }

        residualsKnown = true;
    }

    protected void updateSumOfSquaredResiduals() {
        double delta = 0.0;

        int oldFlag = residualFlag;
        residualFlag = 1 - residualFlag;

        for (int i = 0; i < locationCount; i++) {
            if (locationUpdated[i]) {
                // if location i is updated, calculate the residuals to all js
                // also sum the change in sum residual
                for (int j = 0; j < locationCount; j++) {
                    double distance = calculateDistance(locations[i], locations[j]);
                    double residual = distance - observations[i][j];
                    double squaredResidual = residual * residual;

                    delta += squaredResiduals[oldFlag][i][j] - squaredResidual;

                    squaredResiduals[residualFlag][i][j] = squaredResidual;
                    squaredResiduals[residualFlag][j][i] = squaredResidual;
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

    private double[][] observations;
    private double[][] locations;

    private boolean[] locationUpdated;

    private boolean residualsKnown = false;
    private int residualFlag = 0;
    private int storedResidualFlag;

    private boolean sumOfSquaredResidualsKnown = false;
    private double[][][] squaredResiduals;
    private double sumOfSquaredResiduals;
    private double storedSumOfSquaredResiduals;

    private boolean truncationsKnown = false;
    private double truncationSum;
    private double storedTruncationSum;
    private double[] truncations;
    private double[] storedTruncations;

}

/*
 * MultidimensionalScalingLikelihood.java
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

package dr.evomodel.antigenic;

import dr.inference.model.*;
import dr.math.distributions.NormalDistribution;
import dr.util.DataTable;
import dr.xml.*;

import java.io.FileReader;
import java.io.IOException;

/**
 * @author Andrew Rambaut
 * @author Marc Suchard
 * @version $Id$
 */
public class MultidimensionalScalingLikelihood extends AbstractModelLikelihood {


    public enum ObservationType {
        POINT,
        UPPER_BOUND,
        LOWER_BOUND,
        MISSING
    }

    public final static String MULTIDIMENSIONAL_SCALING_LIKELIHOOD = "multidimensionalScalingLikelihood";

    public MultidimensionalScalingLikelihood(String name) {

        super(name);
    }

    /**
     * A simple constructor for a fully specified symmetrical data matrix
     * @param mdsDimension
     * @param mdsPrecision
     * @param locationsParameter
     * @param dataTable
     */
    public MultidimensionalScalingLikelihood(
            int mdsDimension,
            boolean isLeftTruncated,
            Parameter mdsPrecision,
            MatrixParameter locationsParameter,
            DataTable<double[]> dataTable) {

        super(MULTIDIMENSIONAL_SCALING_LIKELIHOOD);

        // construct a compact data table
        String[] rowLabels = dataTable.getRowLabels();
        String[] columnLabels = dataTable.getRowLabels();

        int rowCount = dataTable.getRowCount();
        int observationCount = ((rowCount - 1) * rowCount) / 2;
        double[] observations = new double[observationCount];
        ObservationType[] observationTypes = new ObservationType[observationCount];
//        int[] distanceIndices = new int[observationCount];
        int[] rowLocationIndices = new int[observationCount];
        int[] columnLocationIndices = new int[observationCount];

        int u = 0;
        for (int i = 0; i < rowCount; i++) {

            double[] dataRow = dataTable.getRow(i);

            for (int j = i + 1; j < rowCount; j++) {
                observations[u] = dataRow[j];
                observationTypes[u] = ObservationType.POINT;
//                distanceIndices[u] = u;
                rowLocationIndices[u] = i;
                columnLocationIndices[u] = j;
                u++;
            }

        }

        initialize(mdsDimension, isLeftTruncated, mdsPrecision, locationsParameter, rowLabels, observations, observationTypes, rowLocationIndices, columnLocationIndices);
    }

    protected void initialize(
            final int mdsDimension,
            final boolean isLeftTruncated,
            final Parameter mdsPrecision,
            final MatrixParameter locationsParameter,
            final String[] locationLabels,
            final double[] observations,
            final ObservationType[] observationTypes,
//            final int[] distanceIndices,
            final int[] rowLocationIndices,
            final int[] columnLocationIndices) {

        this.mdsDimension = mdsDimension;

        locationCount = locationLabels.length;

        // upper triangular matrix
        this.distanceCount = (locationCount * (locationCount - 1)) / 2;

        this.locationLabels = locationLabels;

        this.observations = observations;
        this.observationTypes = observationTypes;
//        this.distanceIndices = distanceIndices;
        this.rowLocationIndices = rowLocationIndices;
        this.columnLocationIndices = columnLocationIndices;

        this.observationCount = observations.length;
        this.upperThresholdCount = 0;
        this.lowerThresholdCount = 0;

        for (ObservationType type : observationTypes) {
            upperThresholdCount += (type == ObservationType.UPPER_BOUND ? 1 : 0);
            lowerThresholdCount += (type == ObservationType.LOWER_BOUND ? 1 : 0);
        }

        thresholdCount = upperThresholdCount + lowerThresholdCount;
        pointObservationCount = observationCount - thresholdCount;

        upperThresholdIndices = new int[upperThresholdCount];
        lowerThresholdIndices = new int[lowerThresholdCount];
        pointObservationIndices = new int[pointObservationCount];

        int ut = 0;
        int lt = 0;
        int po = 0;
        for (int i = 0; i < observationCount; i++) {
            switch (observationTypes[i]) {
                case POINT:
                    pointObservationIndices[po] = i;
                    po++;
                    break;
                case UPPER_BOUND:
                    upperThresholdIndices[ut] = i;
                    ut++;
                    break;
                case LOWER_BOUND:
                    lowerThresholdIndices[lt] = i;
                    lt++;
                    break;
            }
        }

        this.locationsParameter = locationsParameter;
        setupLocationsParameter(this.locationsParameter);
        addVariable(locationsParameter);
        locationUpdated = new boolean[locationsParameter.getParameterCount()];

        // a cache of row to column distances (column indices given by array above).
        distances = new double[distanceCount];
        storedDistances = new double[distanceCount];
        distanceUpdated = new boolean[distanceCount];

        // a cache of individual truncations
        truncations = new double[distanceCount];
        storedTruncations = new double[distanceCount];

        // a cache of threshold calcs
        thresholds = new double[thresholdCount];
        storedThresholds = new double[thresholdCount];

        this.mdsPrecisionParameter = mdsPrecision;
        addVariable(mdsPrecision);

        this.isLeftTruncated = isLeftTruncated; // Re-normalize likelihood for strictly positive distances

        // make sure everything is calculated on first evaluation
        makeDirty();

        addStatistic(new Distances());
    }

    protected void setupLocationsParameter(MatrixParameter locationsParameter) {
        if (locationsParameter.getColumnDimension() > 0){
            boolean proceed = true;
            //TODO: allow for a generic tipTrait/location specification with mismatches: tipTraits which are not represented by the pairwise distances (tipTraits > locations) and more locations based on pairwise distances than represented by the tipTraits
            if (locationsParameter.getColumnDimension() != locationCount){
                System.err.println("locationsParameter column dimension ("+locationsParameter.getColumnDimension()+") is not equal to the locationCount ("+locationCount+")");
                proceed = false;
            }
            if (locationsParameter.getRowDimension() != mdsDimension){
                System.err.println("locationsParameter row dimension ("+locationsParameter.getRowDimension()+") is not equal to the mdsDimension ("+mdsDimension+")");
                proceed = false;
            }
            if (!proceed) {
                System.exit(-1);
            }
        } else{
            locationsParameter.setColumnDimension(mdsDimension);
            locationsParameter.setRowDimension(locationCount);
        }
        for (int i = 0; i < locationLabels.length; i++) {
            locationsParameter.getParameter(i).setId(locationLabels[i]);
        }

        for (int i = 0; i < locationsParameter.getParameterCount(); ++i) {
            Parameter param = locationsParameter.getParameter(i);
            try {
                if (param.getBounds() != null) {
                    // Do nothing
                }
            } catch (NullPointerException exception) {
                param.addBounds(new Parameter.DefaultBounds(
                        Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, param.getDimension()));
            }
        }
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Variable.ChangeType type) {
        // TODO Flag which cachedDistances or mdsPrecision need to be updated

        if (variable == locationsParameter) {
            int location = index / mdsDimension;

            locationUpdated[location] = true;
            distancesKnown = false;
            residualsKnown = false;
            thresholdsKnown = false;
            truncationsKnown = false;
        } else if (variable == mdsPrecisionParameter) {
            for (int i = 0; i < distanceUpdated.length; i++) {
                distanceUpdated[i] = true;
            }
            residualsKnown = false;
            thresholdsKnown = false;
            truncationsKnown = false;
        } else {
            // could be a derived class's parameter
//            throw new IllegalArgumentException("Unknown parameter");
        }
        likelihoodKnown = false;
    }

    @Override
    protected void storeState() {
        System.arraycopy(distances, 0, storedDistances, 0, distances.length);
        System.arraycopy(truncations, 0, storedTruncations, 0, truncations.length);
        System.arraycopy(thresholds, 0, storedThresholds, 0, thresholds.length);

        storedLogLikelihood = logLikelihood;
        storedTruncationSum = truncationSum;
        storedThresholdSum = thresholdSum;
        storedSumOfSquaredResiduals = sumOfSquaredResiduals;
    }

    @Override
    protected void restoreState() {
        double[] tmp = storedDistances;
        storedDistances = distances;
        distances = tmp;
        distancesKnown = true;

        tmp = storedTruncations;
        storedTruncations = truncations;
        truncations = tmp;

        tmp = storedThresholds;
        storedThresholds = thresholds;
        thresholds = tmp;

        logLikelihood = storedLogLikelihood;
        likelihoodKnown = true;

        truncationSum = storedTruncationSum;
        truncationsKnown = true;

        thresholdSum = storedThresholdSum;
        thresholdsKnown = true;

        sumOfSquaredResiduals = storedSumOfSquaredResiduals;
        residualsKnown = true;
    }

    @Override
    protected void acceptState() {
        // do nothing
    }

    public void makeDirty() {
        distancesKnown = false;
        likelihoodKnown = false;
        residualsKnown = false;
        truncationsKnown = false;
        thresholdsKnown = false;

        for (int i = 0; i < locationUpdated.length; i++) {
            locationUpdated[i] = true;
        }
        for (int i = 0; i < distanceUpdated.length; i++) {
            distanceUpdated[i] = true;
        }
    }

    public Model getModel() {
        return this;
    }

    public double getLogLikelihood() {
        if (!likelihoodKnown) {
            if (!distancesKnown) {
                calculateDistances();
                residualsKnown = false;
            }

            logLikelihood = computeLogLikelihood();

            for (int i = 0; i < locationUpdated.length; i++) {
                locationUpdated[i] = false;
            }
            for (int i = 0; i < distanceUpdated.length; i++) {
                distanceUpdated[i] = false;
            }
        }

        return logLikelihood;
    }

    // This function can be overwritten to implement other sampling densities, i.e. discrete ranks
    protected double computeLogLikelihood() {

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

        return logLikelihood;
    }

    protected double calculateThresholdObservations(double precision) {
        double sum = 0.0;
        double sd = 1.0 / Math.sqrt(precision);
        int j = 0;
        for (int i = 0; i < upperThresholdCount; i++) {
            int observationIndex = upperThresholdIndices[i];
            int dist = getDistanceIndexForObservation(observationIndex);
            if (dist != -1) {
                if (distanceUpdated[dist]) {
//                double cdf = NormalDistribution.cdf(observations[observationIndex], distances[dist], sd, false);
//                double tail = 1.0 - cdf;
                    // using special tail function of NormalDistribution (see main() in NormalDistribution for test)
                    double tail = NormalDistribution.tailCDF(observations[observationIndex], distances[dist], sd);
                    thresholds[j] = Math.log(tail);
                }
            } else {
                // -1 denotes a distance to self (i.e., 0)
                double tail = NormalDistribution.tailCDF(observations[observationIndex], 0.0, sd);
                thresholds[j] = Math.log(tail);
            }

            if (Double.isInfinite(thresholds[j])) {
                System.out.println("Error calculation threshold probability");
            }

            sum += thresholds[j];
            j++;
        }
        for (int i = 0; i < lowerThresholdCount; i++) {
            int observationIndex = lowerThresholdIndices[i];
            int dist = getDistanceIndexForObservation(observationIndex);

           if (dist != -1) {
                if (distanceUpdated[dist]) {
                    thresholds[j] = NormalDistribution.cdf(observations[observationIndex], distances[dist], sd, true);
                }
            } else {
                // -1 denotes a distance to self (i.e., 0)
                thresholds[j] = NormalDistribution.cdf(observations[observationIndex], 0.0, sd, true);
            }

            if (Double.isInfinite(thresholds[j])) {
                System.out.println("Error calculation threshold probability");
            }

            sum += thresholds[j];
            j++;
        }

        thresholdsKnown = true;

        return sum;
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

    private int getDistanceIndexForObservation(final int observationIndex) {
        int row = getLocationIndex(rowLocationIndices[observationIndex]);
        int col = getLocationIndex(columnLocationIndices[observationIndex]);

        if (row == col) {
            return -1;
        }

        // make sure row index is less than column index
        if (row > col) {
            int tmp = row;
            row = col;
            col = tmp;
        }

        // get the index of [row, col] in an unrolled upper triangular matrix
        int index = 0;
        for (int x = 0; x < row; x++) {
            index += locationCount - x - 1;
        }
        index += col - row - 1;

        return index;
    }

    /**
     *     Returns the index of the  given by index
     *     Overriding this allows the clustering of locations
     * @param index
     * @return
     */
    protected int getLocationIndex(int index) {
        return index;
    }

    public String[] getLocationLabels() {
        return locationLabels;
    }

    protected double calculateDistance(Parameter X, Parameter Y) {
        double sum = 0.0;
        for (int i = 0; i < mdsDimension; i++) {
            double difference = X.getParameterValue(i) - Y.getParameterValue(i);
            sum += difference * difference;
        }
        return Math.sqrt(sum);
    }

    public int getMDSDimension() {
        return mdsDimension;
    }

    public int getLocationCount() {
        return locationCount;
    }

    public MatrixParameter getLocationsParameter() {
        return locationsParameter;
    }

    public class Distances extends Statistic.Abstract {

        public Distances() {
            super("distances");
        }

        public int getDimension() {
            return distanceCount;
        }

        public double getStatisticValue(int i) {
            if (!distancesKnown) {
                calculateDistances();
            }
            return distances[i];
        }

    }

    // **************************************************************
    // XMLObjectParser
    // **************************************************************

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {
        public final static String FILE_NAME = "fileName";

        public final static String TIP_TRAIT = "tipTrait";
        public final static String LOCATIONS = "locations";
        public static final String MDS_DIMENSION = "mdsDimension";
        public static final String MDS_PRECISION = "mdsPrecision";

        public String getParserName() {
            return MULTIDIMENSIONAL_SCALING_LIKELIHOOD;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            String fileName = xo.getStringAttribute(FILE_NAME);
            DataTable<double[]> distanceTable;
            try {
                distanceTable = DataTable.Double.parse(new FileReader(fileName));
            } catch (IOException e) {
                throw new XMLParseException("Unable to read assay data from file: " + e.getMessage());
            }

            if (distanceTable.getRowCount() != distanceTable.getColumnCount()) {
                throw new XMLParseException("Data table is not symmetrical.");
            }

            int mdsDimension = xo.getIntegerAttribute(MDS_DIMENSION);

            MatrixParameter locationsParameter = (MatrixParameter) xo.getElementFirstChild(LOCATIONS);

            Parameter mdsPrecision = (Parameter) xo.getElementFirstChild(MDS_PRECISION);

            boolean isLeftTruncated = false;

            return new MultidimensionalScalingLikelihood(mdsDimension, isLeftTruncated, mdsPrecision, locationsParameter, distanceTable);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "Provides the likelihood of pairwise distance given vectors of coordinates" +
                    "for points according to the multidimensional scaling scheme of XXX & Rafferty (xxxx).";
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newStringRule(FILE_NAME, false, "The name of the file containing the assay table"),
                AttributeRule.newIntegerRule(MDS_DIMENSION, false, "The dimension of the space for MDS"),
                new ElementRule(LOCATIONS, MatrixParameter.class),
                new ElementRule(MDS_PRECISION, Parameter.class)
        };

        public Class getReturnType() {
            return MultidimensionalScalingLikelihood.class;
        }
    };

    private int distanceCount;
    private int observationCount;
    private int upperThresholdCount;
    private int lowerThresholdCount;
    private int pointObservationCount;
    private int thresholdCount;


    private String[] locationLabels;

    private int locationCount;

    private double[] observations;
    private ObservationType[] observationTypes;
    //    protected int[] distanceIndices;
    private int[] rowLocationIndices;
    private int[] columnLocationIndices;
    private int[] upperThresholdIndices;
    private int[] lowerThresholdIndices;
    private int[] pointObservationIndices;

    private MatrixParameter locationsParameter;
    private Parameter mdsPrecisionParameter;

    private boolean likelihoodKnown = false;
    private double logLikelihood;
    private double storedLogLikelihood;

    protected boolean distancesKnown = false;
    private double sumOfSquaredResiduals;
    private double storedSumOfSquaredResiduals;
    private double[] distances;
    private double[] storedDistances;

    protected boolean[] locationUpdated;
    protected boolean[] distanceUpdated;

    protected boolean residualsKnown = false;

    protected boolean truncationsKnown = false;
    private double truncationSum;
    private double storedTruncationSum;
    private double[] truncations;
    private double[] storedTruncations;

    protected boolean thresholdsKnown = false;
    private double thresholdSum;
    private double storedThresholdSum;
    private double[] thresholds;
    private double[] storedThresholds;

    private boolean isLeftTruncated;
    private int mdsDimension;
}

package dr.evomodel.continuous;

import dr.inference.model.*;
import dr.math.distributions.NormalDistribution;
import dr.util.DataTable;
import dr.xml.*;

import java.io.FileReader;
import java.io.IOException;
import java.util.*;

/**
 * @author Andrew Rambaut
 * @author Marc Suchard
 * @version $Id$
 */
public class MultidimensionalScalingLikelihood extends AbstractModelLikelihood {


    public enum ObservationType {
        POINT,
        UPPER_BOUND,
        LOWER_BOUND
    }

    public final static String MULTIDIMENSIONAL_SCALING_LIKELIHOOD = "multidimensionalScalingLikelihood";

    public MultidimensionalScalingLikelihood(String name) {

        super(name);
    }

    /**
     * A simple constructor for a fully specified symmetrical data matrix
     * @param mdsDimension
     * @param mdsPrecision
     * @param tipTraitParameter
     * @param locationsParameter
     * @param dataTable
     */
    public MultidimensionalScalingLikelihood(
            int mdsDimension,
            Parameter mdsPrecision,
            CompoundParameter tipTraitParameter,
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
        int[] distanceIndices = new int[observationCount];
        int[] rowIndices = new int[observationCount];
        int[] columnIndices = new int[observationCount];

        int u = 0;
        for (int i = 0; i < rowCount; i++) {

            double[] dataRow = dataTable.getRow(i);

            for (int j = i + 1; j < rowCount; j++) {
                observations[u] = dataRow[j];
                observationTypes[u] = ObservationType.POINT;
                distanceIndices[u] = u;
                rowIndices[u] = i;
                columnIndices[u] = j;
                u++;
            }

        }

        initialize(mdsDimension, mdsPrecision, tipTraitParameter, locationsParameter, null, rowLabels, columnLabels, observations, observationTypes, distanceIndices, rowIndices, columnIndices);
    }

    protected void initialize(
            final int mdsDimension,
            final Parameter mdsPrecision,
            final CompoundParameter tipTraitParameter,
            final MatrixParameter rowLocationsParameter,
            final MatrixParameter columnLocationsParameter,
            final String[] rowLabels,
            final String[] columnLabels,
            final double[] observations,
            final ObservationType[] observationTypes,
            final int[] distanceIndices,
            final int[] rowIndices,
            final int[] columnIndices) {

        rowCount = rowLabels.length;
        columnCount = columnLabels.length;

        Map<String, Integer> tipNameMap = null;
        if (tipTraitParameter != null) {
            tipCount = tipTraitParameter.getNumberOfParameters();

            assert(rowCount == tipCount);

            //  the row -> tip map
            tipIndices = new int[tipCount];

            tipNameMap = new HashMap<String, Integer>();
            for (int i = 0; i < tipCount; i++) {
                String label = tipTraitParameter.getParameter(i).getParameterName();
                if (label.endsWith(".antigenic")) {
                    label = label.substring(0, label.indexOf(".antigenic"));
                }
                for (String rowName : rowLabels) {
                    if (label.toUpperCase().startsWith(rowName.toUpperCase())) {
                        break;
                    }
                }

                tipNameMap.put(label, i);

                tipIndices[i] = -1;
            }
        } else {
            tipIndices = null;
            tipCount = 0;
        }

        this.observations = observations;
        this.observationTypes = observationTypes;
        this.distanceIndices = distanceIndices;
        this.rowIndices = rowIndices;
        this.columnIndices = columnIndices;

        for (int i = 0; i < rowCount; i++) {

            if (tipIndices != null) {
                // if the row is in the tree then add a entry to map tip to row
                Integer tipIndex = tipNameMap.get(rowLabels[i]);
                if (tipIndex != null) {
                    tipIndices[tipIndex] = i;
                    rowIndices[i] = tipIndex;
                } else {
                    System.err.println("Tip, " + rowLabels[i] + ", not found in tree");
                }
            }
        }

        this.distancesCount = rowIndices.length;
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

        if (tipIndices != null) {
            for (int i = 0; i < tipCount; i++) {
                if (tipIndices[i] == -1) {
                    String label = tipTraitParameter.getParameter(i).getParameterName();
                    System.err.println("Tree tip, " + label + ", not found in data table");
                }
            }
        }

        // add tipTraitParameter to enable store / restore
        this.tipTraitParameter = tipTraitParameter;
        if (tipTraitParameter != null) {
            addVariable(tipTraitParameter);
        }

        this.rowLocationsParameter = rowLocationsParameter;
        rowLocationsParameter.setColumnDimension(mdsDimension);
        rowLocationsParameter.setRowDimension(rowCount);
        addVariable(rowLocationsParameter);
        rowLocationUpdated = new boolean[rowCount];

        if (columnLocationsParameter != null) {
            this.columnLocationsParameter = columnLocationsParameter;
            columnLocationsParameter.setColumnDimension(mdsDimension);
            columnLocationsParameter.setRowDimension(columnCount);
            addVariable(columnLocationsParameter);
            columnLocationUpdated = new boolean[columnCount];
        } else {
            this.columnLocationsParameter = rowLocationsParameter;
            columnLocationUpdated = rowLocationUpdated;
        }

        // a cache of row to column distances (column indices given by array above).
        distances = new double[distancesCount];
        storedDistances = new double[distancesCount];
        distanceUpdate = new boolean[distancesCount];

        // a cache of individual truncations
        truncations = new double[distancesCount];
        storedTruncations = new double[distancesCount];

        // a cache of threshold calcs
        thresholds = new double[thresholdCount];
        storedThresholds = new double[thresholdCount];

        this.mdsDimension = mdsDimension;

        this.mdsPrecisionParameter = mdsPrecision;
        addVariable(mdsPrecision);

        this.isLeftTruncated = true; // Re-normalize likelihood for strictly positive distances

    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Variable.ChangeType type) {
        // TODO Flag which cachedDistances or mdsPrecision need to be updated

        if (variable == rowLocationsParameter) {
            int rowIndex = index / mdsDimension;
            int dim = index % mdsDimension;

            if (tipTraitParameter != null) {
                if (tipIndices[rowIndex] != -1) {
                    double value = rowLocationsParameter.getParameterValue(index);
                    tipTraitParameter.setParameterValue((rowIndex * mdsDimension) + dim, value);
                }
            }

            rowLocationUpdated[index / mdsDimension] = true;
            distancesKnown = false;
            thresholdsKnown = false;
            truncationKnown = false;
        } else if (variable == columnLocationsParameter) {
            columnLocationUpdated[index / mdsDimension] = true;
            distancesKnown = false;
            thresholdsKnown = false;
            truncationKnown = false;
        } else if (variable == mdsPrecisionParameter) {
            for (int i = 0; i < distanceUpdate.length; i++) {
                distanceUpdate[i] = true;
            }
            thresholdsKnown = false;
            truncationKnown = false;
        } else if (variable == tipTraitParameter) {
//            throw new IllegalArgumentException("Only MultidimensionalScalingLikelihood should be changing the tipTraitParameter");
        } else {
            throw new IllegalArgumentException("Unknown parameter");
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
        truncationKnown = true;

        thresholdSum = storedThresholdSum;
        thresholdsKnown = true;

        sumOfSquaredResiduals = storedSumOfSquaredResiduals;
    }

    @Override
    protected void acceptState() {
        // do nothing
    }

    public void makeDirty() {
        distancesKnown = false;
        likelihoodKnown = false;
        truncationKnown = false;
        thresholdsKnown = false;

        for (int i = 0; i < rowLocationUpdated.length; i++) {
            rowLocationUpdated[i] = true;
        }
        for (int i = 0; i < columnLocationUpdated.length; i++) {
            columnLocationUpdated[i] = true;
        }
        for (int i = 0; i < distanceUpdate.length; i++) {
            distanceUpdate[i] = true;
        }
    }

    public Model getModel() {
        return this;
    }

    public double getLogLikelihood() {
        if (!likelihoodKnown) {
            if (!distancesKnown) {
                calculateDistances();
                sumOfSquaredResiduals = calculateSumOfSquaredResiduals();
                distancesKnown = true;
            }

            logLikelihood = computeLogLikelihood();
            likelihoodKnown = true;
        }

        for (int i = 0; i < rowLocationUpdated.length; i++) {
            rowLocationUpdated[i] = false;
        }
        for (int i = 0; i < columnLocationUpdated.length; i++) {
            columnLocationUpdated[i] = false;
        }
        for (int i = 0; i < distanceUpdate.length; i++) {
            distanceUpdate[i] = false;
        }

        return logLikelihood;
    }

    // This function can be overwritten to implement other sampling densities, i.e. discrete ranks
    protected double computeLogLikelihood() {

        double precision = mdsPrecisionParameter.getParameterValue(0);

        // totalNonMissingCount should be totalObservedCount (not > or < threshold)
        double logLikelihood = (pointObservationCount / 2) * Math.log(precision) - 0.5 * precision * sumOfSquaredResiduals;

        if (thresholdCount > 0) {
            if (!thresholdsKnown) {
                thresholdSum = calculateThresholdObservations(precision);
                thresholdsKnown = true;
            }
            logLikelihood += thresholdSum;
        }

        if (isLeftTruncated) {
            if (!truncationKnown) {
                truncationSum = calculateTruncation(precision);
                truncationKnown = true;
            }
            logLikelihood -= truncationSum;
        }

        return logLikelihood;
    }

    private double calculateThresholdObservations(double precision) {
        double sum = 0.0;
        double sd = 1.0 / Math.sqrt(precision);
        int j = 0;
        for (int i = 0; i < upperThresholdCount; i++) {
            int observationIndex = upperThresholdIndices[i];
            int dist = distanceIndices[observationIndex];
            if (distanceUpdate[dist]) {
                thresholds[j] = Math.log(1.0 - NormalDistribution.cdf(observations[observationIndex], distances[dist], sd));
            }
            sum += thresholds[j];
            j++;
        }
        for (int i = 0; i < lowerThresholdCount; i++) {
            int observationIndex = upperThresholdIndices[i];
            int dist = distanceIndices[observationIndex];
            if (distanceUpdate[dist]) {
                thresholds[j] = Math.log(NormalDistribution.cdf(observations[observationIndex], distances[dist], sd));
            }
            sum += thresholds[j];
            j++;
        }

        return sum;
    }

    private double calculateTruncation(double precision) {
        double sum = 0.0;
        double sd = 1.0 / Math.sqrt(precision);
        for (int i = 0; i < observationCount; i++) {
            if (distanceUpdate[i]) {
                truncations[i] = Math.log(NormalDistribution.cdf(distances[distanceIndices[i]], 0.0, sd));
            }
            sum += truncations[i];
        }
        return sum;
    }

    private double calculateSumOfSquaredResiduals() {
        double sum = 0.0;
        for (int i = 0; i < observationCount; i++) {
            if (observationTypes[i] == ObservationType.POINT) {
                // Only increment sum if dataTable[i][j] is observed (not > or < threshold)
                double residual = distances[distanceIndices[i]] - observations[i];
                sum += residual * residual;
            }
        }
        return sum;
    }

    private void calculateDistances() {
        for (int i = 0; i < distancesCount; i++) {
            if (rowLocationUpdated[rowIndices[i]] || columnLocationUpdated[columnIndices[i]]) {
                distances[i] = calculateDistance(
                        rowLocationsParameter.getParameter(rowIndices[i]),
                        columnLocationsParameter.getParameter(columnIndices[i]));
                distanceUpdate[i] = true;
            }
        }
    }

    private double calculateDistance(Parameter X, Parameter Y) {
        double sum = 0.0;
        for (int i = 0; i < mdsDimension; i++) {
            double difference = X.getParameterValue(i) - Y.getParameterValue(i);
            sum += difference * difference;
        }
        return Math.sqrt(sum);
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
                throw new XMLParseException("Unable to read assay data from file, " + fileName);
            }

            if (distanceTable.getRowCount() != distanceTable.getColumnCount()) {
                throw new XMLParseException("Data table is not symmetrical.");
            }

            int mdsDimension = xo.getIntegerAttribute(MDS_DIMENSION);

            // This parameter needs to be linked to the one in the IntegratedMultivariateTreeLikelihood (I suggest that the parameter is created
            // here and then a reference passed to IMTL - which optionally takes the parameter of tip trait values, in which case it listens and
            // updates accordingly.
            CompoundParameter tipTraitParameter = null;
            if (xo.hasChildNamed(TIP_TRAIT)) {
                tipTraitParameter = (CompoundParameter) xo.getElementFirstChild(TIP_TRAIT);
            }

            MatrixParameter locationsParameter = (MatrixParameter) xo.getElementFirstChild(LOCATIONS);

            Parameter mdsPrecision = (Parameter) xo.getElementFirstChild(MDS_PRECISION);

            if (tipTraitParameter != null) {
                if (distanceTable.getRowCount() != tipTraitParameter.getNumberOfParameters()) {
                    throw new XMLParseException("Tree has different number of tips than the number of rows in the distance matrix");
                }
            }

            return new MultidimensionalScalingLikelihood(mdsDimension, mdsPrecision, tipTraitParameter, locationsParameter, distanceTable);
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
                new ElementRule(TIP_TRAIT, CompoundParameter.class, "The parameter of tip locations from the tree", true),
                new ElementRule(LOCATIONS, MatrixParameter.class),
                new ElementRule(MDS_PRECISION, Parameter.class)
        };

        public Class getReturnType() {
            return MultidimensionalScalingLikelihood.class;
        }
    };

    private int distancesCount;
    private int observationCount;
    private int upperThresholdCount;
    private int lowerThresholdCount;
    private int pointObservationCount;
    private int thresholdCount;

    private int tipCount;
    private int rowCount;
    private int columnCount;

    private double[] observations;
    private ObservationType[] observationTypes;
    private int[] distanceIndices;
    private int[] rowIndices;
    private int[] columnIndices;
    private int[] tipIndices;
    private int[] upperThresholdIndices;
    private int[] lowerThresholdIndices;
    private int[] pointObservationIndices;

    private CompoundParameter tipTraitParameter;
    private MatrixParameter rowLocationsParameter;
    private MatrixParameter columnLocationsParameter;
    private Parameter mdsPrecisionParameter;

    private boolean likelihoodKnown = false;
    private double logLikelihood;
    private double storedLogLikelihood;

    private boolean distancesKnown = false;
    private double sumOfSquaredResiduals;
    private double storedSumOfSquaredResiduals;
    private double[] distances;
    private double[] storedDistances;

    private boolean[] rowLocationUpdated;
    private boolean[] columnLocationUpdated;
    private boolean[] distanceUpdate;

    private boolean truncationKnown = false;
    private double truncationSum;
    private double storedTruncationSum;
    private double[] truncations;
    private double[] storedTruncations;

    private boolean thresholdsKnown = false;
    private double thresholdSum;
    private double storedThresholdSum;
    private double[] thresholds;
    private double[] storedThresholds;

    private boolean isLeftTruncated;
    private int mdsDimension;
}

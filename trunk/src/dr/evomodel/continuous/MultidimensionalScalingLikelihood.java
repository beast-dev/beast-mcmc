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

    public final static String MULTIDIMENSIONAL_SCALING_LIKELIHOOD = "multidimensionalScalingLikelihood";

    public MultidimensionalScalingLikelihood(String name) {

        super(name);
    }

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
        double[][][] compactDataMatrix = new double[rowCount][][];

        int[][] nonMissingIndices = new int[rowCount][];

        for (int i = 0; i < rowCount; i++) {

            double[] dataRow = dataTable.getRow(i);

            int nonMissingCount = 0;
            for (int j = 0; j < rowCount; j++) {
                if (!Double.isNaN(dataRow[j]) && dataRow[j] > 0) {
                    nonMissingCount ++;
                }
            }

            compactDataMatrix[i] = new double[nonMissingCount][1];
            nonMissingIndices[i] = new int[nonMissingCount];

            int k = 0;
            for (int j = 0; j < rowCount; j++) {
                if (!Double.isNaN(dataRow[j]) && dataRow[j] > 0) {
                    compactDataMatrix[i][k][0] = dataRow[j];
                    nonMissingIndices[i][k] = j;
                    k ++;
                }
            }
            totalNonMissingCount += nonMissingCount;
        }

        initialize(mdsDimension, mdsPrecision, tipTraitParameter, locationsParameter, null, rowLabels, columnLabels, compactDataMatrix, nonMissingIndices, null);
    }

    protected void initialize(
            final int mdsDimension,
            final Parameter mdsPrecision,
            final CompoundParameter tipTraitParameter,
            final MatrixParameter rowLocationsParameter,
            final MatrixParameter columnLocationsParameter,
            final String[] rowLabels,
            final String[] columnLabels,
            final double[][][] compactDataMatrix,
            final int[][] nonMissingIndices,
            boolean[][][] isThreshold) {

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
                for (String rowName : rowLabels) {
                    if (label.toUpperCase().equals(rowName.toUpperCase())) {
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

        this.compactDataMatrix = compactDataMatrix;
        this.nonMissingIndices = nonMissingIndices;

        if (tipIndices != null) {
            rowIndices = new int[rowCount];
        } else {
            rowIndices = null;
        }

        this.isThreshold = isThreshold;
        hasThresholds = (isThreshold != null);

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

        this.totalDistancesCount = 0;
        this.totalNonMissingCount = 0;

        for (int i = 0; i < compactDataMatrix.length; i++) {

            totalDistancesCount += compactDataMatrix[i].length;

            for (int j = 0; j < compactDataMatrix[i].length; j++) {
                for (int k = 0; k < isThreshold[i][j].length; k++) {
                    if (!isThreshold[i][j][k]) {
                        // only count as measured data if not a threshold
                        totalNonMissingCount ++;
                    }
                }
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
        distances = new double[totalNonMissingCount];
        storedDistances = new double[totalNonMissingCount];

        distanceUpdate = new boolean[totalNonMissingCount];

        // a cache of individual truncations
        truncations = new double[totalNonMissingCount];
        storedTruncations = new double[totalNonMissingCount];

        this.mdsDimension = mdsDimension;

        this.mdsPrecisionParameter = mdsPrecision;
        addVariable(mdsPrecision);

        this.isLeftTruncated = false; // Re-normalize likelihood for strictly positive distances

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

            makeDirty();
        } else if (variable == columnLocationsParameter) {
            columnLocationUpdated[index / mdsDimension] = true;
            distancesKnown = false;

            makeDirty();
        } else if (variable == mdsPrecisionParameter) {
            for (int i = 0; i < distanceUpdate.length; i++) {
                distanceUpdate[i] = true;
            }
        } else if (variable == tipTraitParameter) {
            throw new IllegalArgumentException("Only MultidimensionalScalingLikelihood should be changing the tipTraitParameter");
        } else {
            throw new IllegalArgumentException("Unknown parameter");
        }

        truncationKnown = false;
        likelihoodKnown = false;
    }

    @Override
    protected void storeState() {
        System.arraycopy(distances, 0, storedDistances, 0, distances.length);
        System.arraycopy(truncations, 0, storedTruncations, 0, truncations.length);

        storedLogLikelihood = logLikelihood;
        storedTruncationSum = truncationSum;
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

        logLikelihood = storedLogLikelihood;
        likelihoodKnown = true;

        truncationSum = storedTruncationSum;
        truncationKnown = true;

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
        double logLikelihood = (totalNonMissingCount / 2) * Math.log(precision) - 0.5 * precision * sumOfSquaredResiduals;

        if (hasThresholds) {
            logLikelihood += calculateThresholdedObservations(precision);
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

    private final static double minThresholdValue = 0.0;

    private double calculateThresholdedObservations(double precision) {
        double logProbability = 0.0;
        double sd = 1.0 / Math.sqrt(precision);
        int u = 0;
        int v = 0;
        for (int i = 0; i < compactDataMatrix.length; i++) {
            for (int j = 0; j < compactDataMatrix[i].length; j++) {
                if (distanceUpdate[u]) {
                    for (int k = 0; k < compactDataMatrix[i][j].length; k++) {
                        if (isThreshold[i][j][k]) {
                            // TODO Check: switch minThresholdValue and distances[k] order?

                            thresholds[v] = Math.log(NormalDistribution.cdf(minThresholdValue, distances[k], sd));
                        } else {
                            thresholds[v] = 0.0;
                        }
                    }
                    v++;
                }
                u++;
            }
        }

        // TODO Check: + or - thresholds[k]?        
        for (int k = 0; k < thresholds.length; k++) {
            logProbability += thresholds[k];
        }

        return logProbability;
    }

    private double calculateTruncation(double precision) {
        double sum = 0.0;
        double sd = 1.0 / Math.sqrt(precision);
        int u = 0;
        int v = 0;
        for (int i = 0; i < compactDataMatrix.length; i++) {
            for (int j = 0; j < compactDataMatrix[i].length; j++) {
                if (distanceUpdate[u]) {
                    for (int k = 0; k < compactDataMatrix[i][j].length; k++) {
                        truncations[v] = Math.log(NormalDistribution.cdf(distances[k], 0.0, sd));
                        v++;
                    }
                }
                u++;
            }
        }

        for (int k = 0; k < truncations.length; k++) {
            sum += truncations[k];
        }

        return sum;
    }

    private double calculateSumOfSquaredResiduals() {
        double sum = 0.0;
        int u = 0;
        for (int i = 0; i < compactDataMatrix.length; i++) {
            for (int j = 0; j < compactDataMatrix[i].length; j++) {
                for (int k = 0; k < compactDataMatrix[i][j].length; k++) {
                    if (isThreshold[i][j][k]) {
                        // Only increment sum if dataTable[i][j] is observed (not > or < threshold)
                        double residual = distances[u] - compactDataMatrix[i][j][k];
                        sum += residual * residual;
                    }
                }
                u++;
            }
        }
        return sum;
    }

    private void calculateDistances() {
        int k = 0;
        for (int i = 0; i < compactDataMatrix.length; i++) {
            if (rowLocationUpdated[i]) {
                for (int j = 0; j < compactDataMatrix[i].length; j++) {
                    if (columnLocationUpdated[i]) {
                        distances[k] = calculateDistance(rowLocationsParameter.getParameter(i),
                                columnLocationsParameter.getParameter(nonMissingIndices[i][j]));
                        distanceUpdate[k] = true;
                    }
                    k++;
                }
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

    private double[][][] compactDataMatrix;

    private int tipCount;
    private int rowCount;
    private int columnCount;
    private int totalDistancesCount;
    private int totalNonMissingCount;
    private int[] tipIndices;
    private int[] rowIndices;
    private int[][] nonMissingIndices;

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

    private boolean[][][] isThreshold;
    private double thresholdSum;
    private double storedThresholdSum;
    private double[] thresholds;
    private double[] storedThresholds;

    private boolean isLeftTruncated;
    private int mdsDimension;

    private boolean hasThresholds = false;
}

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

    public MultidimensionalScalingLikelihood(
            int mdsDimension,
            Parameter mdsPrecision,
            CompoundParameter tipTraitParameter,
            MatrixParameter locationsParameter,
            DataTable<double[]> dataTable) {

        super(MULTIDIMENSIONAL_SCALING_LIKELIHOOD);

        this.mdsDimension = mdsDimension;

        String[] rowNames = dataTable.getRowLabels();

        rowCount = dataTable.getRowCount();

        Map<String, Integer> tipNameMap = null;
        if (tipTraitParameter != null) {
            tipCount = tipTraitParameter.getNumberOfParameters();

            assert(rowCount == tipCount);

            //  the row -> tip map
            tipIndices = new int[tipCount];

            tipNameMap = new HashMap<String, Integer>();
            for (int i = 0; i < tipCount; i++) {
                String label = tipTraitParameter.getParameter(i).getParameterName();
                for (String rowName : rowNames) {
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

        locations = new double[rowCount][mdsDimension];

        String[] rowLabels = dataTable.getRowLabels();
        this.dataTable = new double[rowCount][];

        if (tipIndices != null) {
            rowIndices = new int[rowCount];
        } else {
            rowIndices = null;
        }

        int totalNonMissingCount = 0;
        nonMissingIndices = new int[rowCount][];

        for (int i = 0; i < rowCount; i++) {

            double[] dataRow = dataTable.getRow(i);

            if (tipIndices != null) {
                // if the virus is in the tree then add a entry to map tip to virus
                Integer tipIndex = tipNameMap.get(rowLabels[i]);
                if (tipIndex != null) {
                    tipIndices[tipIndex] = i;
                    rowIndices[i] = tipIndex;
                } else {
                    System.err.println("Tip, " + rowLabels[i] + ", not found in tree");
                }
            }

            int nonMissingCount = 0;
            for (int j = 0; j < rowCount; j++) {
                if (!Double.isNaN(dataRow[j]) && dataRow[j] > 0) {
                    nonMissingCount ++;
                }
            }

            this.dataTable[i] = new double[nonMissingCount];
            nonMissingIndices[i] = new int[nonMissingCount];

            int k = 0;
            for (int j = 0; j < rowCount; j++) {
                if (!Double.isNaN(dataRow[j]) && dataRow[j] > 0) {
                    this.dataTable[i][k] = dataRow[j];
                    nonMissingIndices[i][k] = j;
                    k ++;
                }
            }
            totalNonMissingCount += nonMissingCount;
        }

        this.totalNonMissingCount = totalNonMissingCount;

        // a cache of virus to serum distances (serum indices given by array above).
        distances = new double[totalNonMissingCount];
        storedDistances = new double[totalNonMissingCount];

        locationUpdates = new boolean[rowCount];
        distanceUpdate = new boolean[totalNonMissingCount];

        // a cache of individual truncations
        truncations = new double[totalNonMissingCount];
        storedTruncations = new double[totalNonMissingCount];

        if (tipIndices != null) {
            for (int i = 0; i < tipCount; i++) {
                if (tipIndices[i] == -1) {
                    String label = tipTraitParameter.getParameter(i).getParameterName();
                    System.err.println("Tree tip, " + label + ", not found in virus assay table");
                }
            }
        }

        // add tipTraitParameter to enable store / restore
        this.tipTraitParameter = tipTraitParameter;
        if (tipTraitParameter != null) {
            addVariable(tipTraitParameter);
        }

        this.locationsParameter = locationsParameter;
        locationsParameter.setColumnDimension(mdsDimension);
        locationsParameter.setRowDimension(rowCount);
        addVariable(locationsParameter);

        // some random initial locations
//        for (int i = 0; i < virusCount; i++) {
//            locationsParameter.getParameter(i).setId(rowNames[i]);
//            for (int j = 0; j < mdsDimension; j++) {
//                double r = MathUtils.nextGaussian();
//                locationsParameter.getParameter(i).setParameterValue(j, r);
//            }
//        }

        this.mdsPrecisionParameter = mdsPrecision;
        addVariable(mdsPrecision);

        this.isLeftTruncated = false; // Re-normalize likelihood for strictly positive distances

        addStatistic(meanStatistic);
    }

    private double transform(final double value, final double maxValue) {
        // transform to log_2
        double t =  Math.log(maxValue / value) / Math.log(2.0);
        return t;
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Variable.ChangeType type) {
        // TODO Flag which cachedDistances or mdsPrecision need to be updated

        if (variable == locationsParameter) {
            int rowIndex = index / mdsDimension;
            int dim = index % mdsDimension;

            if (tipTraitParameter != null) {
                if (tipIndices[rowIndex] != -1) {
                    double value = locationsParameter.getParameterValue(index);
                    tipTraitParameter.setParameterValue((rowIndex * mdsDimension) + dim, value);
                }
            }

            locationUpdates[index / mdsDimension] = true;
            distancesKnown = false;

            statsKnown = false;

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

        statsKnown = false;
    }

    @Override
    protected void acceptState() {
        // do nothing
    }

    public void makeDirty() {
        distancesKnown = false;
        likelihoodKnown = false;
        truncationKnown = false;

        for (int i = 0; i < locationUpdates.length; i++) {
            locationUpdates[i] = true;
        }

        for (int i = 0; i < distanceUpdate.length; i++) {
            distanceUpdate[i] = true;
        }
    }

    public Model getModel() {
        return this;
    }

    public double getLogLikelihood() {
        makeDirty();
        if (!likelihoodKnown) {
            if (!distancesKnown) {
                calculateDistances();
                sumOfSquaredResiduals = calculateSumOfSquaredResiduals();
                distancesKnown = true;

            }

            logLikelihood = computeLogLikelihood();
            likelihoodKnown = true;
        }

        for (int i = 0; i < locationUpdates.length; i++) {
            locationUpdates[i] = false;
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

        if (hasThresholdedValues) {
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

    private double calculateThresholdedObservations(double precision) {
        double logProbability = 0.0;
        double sd = 1.0 / Math.sqrt(precision);
        int k = 0;
        for (int i = 0; i < dataTable.length; i++) {
            for (int j = 0; j < dataTable[i].length; j++) {
                if (distanceUpdate[k]) {
                    if (isThresholded[k]) {
                        // TODO Check: switch minThresholdValue and distances[k] order?
                        thresholds[k] = Math.log(NormalDistribution.cdf(minThresholdValue, distances[k], sd));
                    } else {
                        thresholds[k] = 0.0;
                    }
                }
                k++;
            }
        }

        // TODO Check: + or - thresholds[k]?        
        for (k = 0; k < thresholds.length; k++) {
            logProbability += thresholds[k];
        }

        return logProbability;
    }

    private double calculateTruncation(double precision) {
        double sum = 0.0;
        double sd = 1.0 / Math.sqrt(precision);
        int k = 0;
        for (int i = 0; i < dataTable.length; i++) {
            for (int j = 0; j < dataTable[i].length; j++) {
                if (distanceUpdate[k]) {
                    truncations[k] = Math.log(NormalDistribution.cdf(distances[k], 0.0, sd));
                }
                k++;
            }
        }

        for ( k = 0; k < truncations.length; k++) {
            sum += truncations[k];
        }

        return sum;
    }

    private double calculateSumOfSquaredResiduals() {
        double sum = 0.0;
        int k = 0;
        for (int i = 0; i < dataTable.length; i++) {
            for (int j = 0; j < dataTable[i].length; j++) {
                // Only increment sum if dataTable[i][j] is observed (not > or < threshold)
                double residual = distances[k] - dataTable[i][j];
                sum += residual * residual;
                k++;
            }
        }
        return sum;
    }

    private void calculateDistances() {
        int k = 0;
        for (int i = 0; i < dataTable.length; i++) {
            for (int j = 0; j < dataTable[i].length; j++) {
                if (locationUpdates[i]) {
                    distances[k] = calculateDistance(locationsParameter.getParameter(i),
                            locationsParameter.getParameter(nonMissingIndices[i][j]));
                    distanceUpdate[k] = true;
                }
                k++;
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

    private final double[][] locations;
    private boolean statsKnown = false;

    private void calculateStats() {
        locationMean = new double[mdsDimension];

        for (int i = 0; i < rowCount; i++) {
            for (int j = 0; j < mdsDimension; j++) {
                locations[i][j] = locationsParameter.getParameter(i).getParameterValue(j);
                locationMean[j] += locations[i][j];
            }
        }
        for (int j = 0; j < mdsDimension; j++) {
            locationMean[j] /= rowCount;
        }

        statsKnown = true;
    }

    private final Statistic meanStatistic = new Statistic.Abstract() {

        public String getStatisticName() {
            return "mean";
        }

        public int getDimension() {
            return mdsDimension;
        }

        public double getStatisticValue(int dim) {
            if (!statsKnown) {
                calculateStats();
            }
            return locationMean[dim];
        }

    };

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

    private final double[][] dataTable;

    private final int tipCount;
    private final int rowCount;
    private final int totalNonMissingCount;
    private final int[] tipIndices;
    private final int[] rowIndices;
    private final int[][] nonMissingIndices;


    private final CompoundParameter tipTraitParameter;
    private final MatrixParameter locationsParameter;
    private final Parameter mdsPrecisionParameter;

    private boolean likelihoodKnown = false;
    private double logLikelihood;
    private double storedLogLikelihood;

    private double[] locationMean;

    private boolean distancesKnown = false;
    private double sumOfSquaredResiduals;
    private double storedSumOfSquaredResiduals;
    private double[] distances;
    private double[] storedDistances;

    private boolean[] locationUpdates;
    private boolean[] distanceUpdate;

    private boolean truncationKnown = false;
    private double truncationSum;
    private double storedTruncationSum;
    private double[] truncations;
    private double[] storedTruncations;

    private boolean[] isThresholded;
    private double thresholdSum;
    private double storedThresholdSum;
    private double[] thresholds;
    private double[] storedThresholds;

    private final boolean isLeftTruncated;
    private final int mdsDimension;

    private boolean hasThresholdedValues = false;
    private double minThresholdValue = 20.0; // TODO Transform
}

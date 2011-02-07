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
public class AntigenicTraitLikelihood extends AbstractModelLikelihood {

    public final static String ANTIGENIC_TRAIT_LIKELIHOOD = "antigenicTraitLikelihood";

    public AntigenicTraitLikelihood(
            int mdsDimension,
            Parameter mdsPrecision,
            CompoundParameter tipTraitParameter,
            MatrixParameter virusLocationsParameter,
            MatrixParameter serumLocationsParameter,
            DataTable<String[]> dataTable,
            final boolean log2Transform,
            final double threshold) {

        super(ANTIGENIC_TRAIT_LIKELIHOOD);

        this.mdsDimension = mdsDimension;

        this.titrationThreshold = threshold;

        String[] virusNames = dataTable.getRowLabels();
        boolean hasAliases = false;
        if (virusNames[0].equals("aliases")) {
            // if the first row is labelled aliases then there are aliases so the
            // viruses start on row 2.
            hasAliases = true;
            virusNames = Arrays.copyOfRange(virusNames, 1, virusNames.length);
        }

//        mdsDimension = virusLocationsParameter.getColumnDimension();

        // the total number of viruses is the number of rows in the table
        virusLocationCount = dataTable.getRowCount() - (hasAliases ? 1 : 0);
        int assayCount = dataTable.getColumnCount();
        tipCount = virusLocationCount;

        String[] serumNames = dataTable.getColumnLabels();
        int[] aliasIndices = aliasIndices = new int[serumNames.length];

        Set<String> aliasSet = null;
        List<String> aliasNames = null;

        if (hasAliases) {
            aliasSet = new HashSet<String>();
            aliasNames = new ArrayList<String>();

            String[] aliases = dataTable.getRow(0);
            for (int i = 0; i < serumNames.length; i++) {
                if (aliasSet.contains(aliases[i])) {
                    aliasIndices[i] = aliasNames.indexOf(aliases[i]);
                } else {
                    aliasSet.add(aliases[i]);
                    aliasNames.add(aliases[i]);
                    aliasIndices[i] = aliasNames.size() - 1;
                }
            }

            // the number of serum locations is the number of aliases
            serumLocationCount = aliasNames.size();
        } else {
            // the number of serum locations is the number of columns
            serumLocationCount = assayCount;

            // one alias for one serum
            for (int i = 0; i < serumNames.length; i++) {
                aliasIndices[i] = i;
            }
        }


        Map<String, Integer> tipNameMap = null;
        if (tipTraitParameter != null) {
            if (tipCount != tipTraitParameter.getNumberOfParameters()) {
                System.err.println("Tree has different number of tips than the number of viruses");
            }

            //  the tip -> virus map
            tipIndices = new int[tipCount];

            tipNameMap = new HashMap<String, Integer>();
            for (int i = 0; i < tipCount; i++) {
                String label = tipTraitParameter.getParameter(i).getParameterName();
                for (String virus : virusNames) {
                    if (label.startsWith(virus)) {
                        label = virus;
                        break;
                    }
                }

                tipNameMap.put(label, i);

                tipIndices[i] = -1;
            }
        } else {
            tipIndices = null;
        }

        // the virus -> tip map
        virusIndices = new int[virusLocationCount];

        locations = new double[virusLocationCount][mdsDimension];
//        locations = new double[serumCount][mdsDimension];

        // a set of vectors for each virus giving serum indices for which assay data is available
        measuredSerumIndices = new int[virusLocationCount][];

        // a compressed (no missing values) set of measured assay values between virus and sera.
        this.assayTable = new double[virusLocationCount][];

        double[] maxSerum = null;

        if (log2Transform) {
            maxSerum = new double[serumLocationCount];
            for (int j = 0; j < assayCount; j++) {
                maxSerum[j] = 0;
                String[] dataColumn = dataTable.getColumn(j);
                for (int i = 0; i < virusLocationCount; i++) {
                    double value = convertString(dataColumn[i]);
                    if (!Double.isNaN(value) && value > 0) {
                        if (value > maxSerum[j]) {
                            maxSerum[j] = value;
                        }
                    }
                }
            }
        }

        int totalMeasurementCount = 0;
        int start = hasAliases ? 1 : 0;

        // Build a sparse matrix of assay values
        for (int i = start; i < virusLocationCount; i++) {
            virusIndices[i] = -1;

            String[] dataRow = dataTable.getRow(i);

            if (tipIndices != null) {
                // if the virus is in the tree then add a entry to map tip to virus
                Integer tipIndex = tipNameMap.get(virusNames[i]);
                if (tipIndex != null) {
                    tipIndices[tipIndex] = i;
                    virusIndices[i] = tipIndex;
                } else {
                    System.err.println("Virus, " + virusNames[i] + ", not found in tree");
                }
            }

            // count the number of measured assay values
            int measuredCount = 0;
            for (int j = 0; j < assayCount; j++) {
                double value = convertString(dataRow[i]);
                if (!Double.isNaN(value) &&value > 0) {
                    measuredCount ++;
                }
            }

            assayTable[i] = new double[measuredCount];
            measuredSerumIndices[i] = new int[measuredCount];

            // now fill in the matrix and keep track of indices
            int k = 0;
            for (int j = 0; j < assayCount; j++) {
                double value = convertString(dataRow[i]);
                if (!Double.isNaN(value) && value > 0) {
                    if (log2Transform) {
                        if (value < titrationThreshold) {
                            this.assayTable[i][k] = Double.POSITIVE_INFINITY;
                        } else {
                            this.assayTable[i][k] = transform(value, maxSerum[j]);
                        }
                    } else {
                        this.assayTable[i][k] = value;
                    }
                    measuredSerumIndices[i][k] = aliasIndices[j];
                    k ++;
                }
            }
            totalMeasurementCount += measuredCount;
        }

        this.totalMeasurementCount = totalMeasurementCount;

        // a cache of virus to serum distances (serum indices given by array above).
        distances = new double[totalMeasurementCount];
        storedDistances = new double[totalMeasurementCount];

        virusLocationUpdates = new boolean[virusLocationCount];
        serumLocationUpdates = new boolean[serumLocationCount];
        distanceUpdate = new boolean[totalMeasurementCount];

        // a cache of individual truncations
        truncations = new double[totalMeasurementCount];
        storedTruncations = new double[totalMeasurementCount];

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

        this.virusLocationsParameter = virusLocationsParameter;
        virusLocationsParameter.setColumnDimension(mdsDimension);
        virusLocationsParameter.setRowDimension(virusLocationCount);
        addVariable(virusLocationsParameter);

        // some random initial locations
//        for (int i = 0; i < virusCount; i++) {
//            virusLocationsParameter.getParameter(i).setId(virusNames[i]);
//            for (int j = 0; j < mdsDimension; j++) {
//                double r = MathUtils.nextGaussian();
//                virusLocationsParameter.getParameter(i).setParameterValue(j, r);
//            }
//        }

        if (serumLocationsParameter != null) {
            this.serumLocationsParameter = serumLocationsParameter;
            serumLocationsParameter.setColumnDimension(mdsDimension);
            serumLocationsParameter.setRowDimension(serumLocationCount);
            addVariable(serumLocationsParameter);

            // some random initial locations
//            for (int i = 0; i < serumCount; i++) {
//                serumLocationsParameter.getParameter(i).setId(serumNames[i]);
//                for (int j = 0; j < mdsDimension; j++) {
//                    double r = MathUtils.nextGaussian();
//                    serumLocationsParameter.getParameter(i).setParameterValue(j, r);
//                }
//            }
        } else {
            this.serumLocationsParameter = virusLocationsParameter;
        }

        this.mdsParameter = mdsPrecision;
        addVariable(mdsPrecision);

        this.isLeftTruncated = false; // Re-normalize likelihood for strictly positive distances

        addStatistic(meanStatistic);
    }

    private double convertString(String value) {
        try {
            return java.lang.Double.valueOf(value);
        } catch (NumberFormatException nfe) {
            return java.lang.Double.NaN;
        }
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

        if (variable == virusLocationsParameter) {
            int virusIndex = index / mdsDimension;
            int dim = index % mdsDimension;

            if (tipTraitParameter != null) {
                if (tipIndices[virusIndex] != -1) {
                    double value = virusLocationsParameter.getParameterValue(index);
                    tipTraitParameter.setParameterValue((virusIndex * mdsDimension) + dim, value);
                }
            }

            virusLocationUpdates[index / mdsDimension] = true;
            distancesKnown = false;

            statsKnown = false;

            makeDirty();
        } else if (variable == serumLocationsParameter) {
            serumLocationUpdates[index / mdsDimension] = true;

            distancesKnown = false;

        } else if (variable == mdsParameter) {
            for (int i = 0; i < distanceUpdate.length; i++) {
                distanceUpdate[i] = true;
            }
        } else if (variable == tipTraitParameter) {
            // do nothing
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

        for (int i = 0; i < virusLocationUpdates.length; i++) {
            virusLocationUpdates[i] = true;
        }

        for (int i = 0; i < serumLocationUpdates.length; i++) {
            serumLocationUpdates[i] = true;
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

        for (int i = 0; i < virusLocationUpdates.length; i++) {
            virusLocationUpdates[i] = false;
        }

        for (int i = 0; i < serumLocationUpdates.length; i++) {
            serumLocationUpdates[i] = false;
        }

        for (int i = 0; i < distanceUpdate.length; i++) {
            distanceUpdate[i] = false;
        }

        return logLikelihood;
    }

    // This function can be overwritten to implement other sampling densities, i.e. discrete ranks
    protected double computeLogLikelihood() {

        double precision = mdsParameter.getParameterValue(0);
        double logLikelihood = (totalMeasurementCount / 2) * Math.log(precision) - 0.5 * precision * sumOfSquaredResiduals;

        if (isLeftTruncated) {
            if (!truncationKnown) {
                truncationSum = calculateTruncation(precision);
                truncationKnown = true;
            }
            logLikelihood -= truncationSum;
        }

        return logLikelihood;
    }

    private double calculateTruncation(double precision) {
        double sum = 0.0;
        double sd = 1.0 / Math.sqrt(precision);
        int k = 0;
        for (int i = 0; i < assayTable.length; i++) {
            for (int j = 0; j < assayTable[i].length; j++) {
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
        for (int i = 0; i < assayTable.length; i++) {
            for (int j = 0; j < assayTable[i].length; j++) {
                double residual = distances[k] - assayTable[i][j];
                sum += residual * residual;
                k++;
            }
        }
        return sum;
    }

    private void calculateDistances() {
        int k = 0;
        for (int i = 0; i < assayTable.length; i++) {
            for (int j = 0; j < assayTable[i].length; j++) {
                if (virusLocationUpdates[i] || serumLocationUpdates[measuredSerumIndices[i][j]]) {
                    distances[k] = calculateDistance(virusLocationsParameter.getParameter(i),
                            serumLocationsParameter.getParameter(measuredSerumIndices[i][j]));
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

        for (int i = 0; i < virusLocationCount; i++) {
            for (int j = 0; j < mdsDimension; j++) {
                locations[i][j] = virusLocationsParameter.getParameter(i).getParameterValue(j);
                locationMean[j] += locations[i][j];
            }
        }
        for (int j = 0; j < mdsDimension; j++) {
            locationMean[j] /= virusLocationCount;
        }

//        for (int i = 0; i < virusCount; i++) {
//            for (int j = 0; j < mdsDimension; j++) {
//                locations[i][j] -= locationMean[j];
//            }
//        }
//
////        for (int i = 0; i < serumCount; i++) {
////            for (int j = 0; j < mdsDimension; j++) {
////                locations[i][j] = serumLocationsParameter.getParameter(i).getParameterValue(j);
////                locationMean[j] += locations[i][j];
////            }
////        }
////        for (int j = 0; j < mdsDimension; j++) {
////            locationMean[j] /= serumCount;
////        }
////
////        for (int i = 0; i < serumCount; i++) {
////            for (int j = 0; j < mdsDimension; j++) {
////                locations[i][j] -= locationMean[j];
////            }
////        }
//
//        RealMatrix data = MatrixUtils.createRealMatrix(locations);
//        // compute the covariance matrix
//        RealMatrix covMatrix = null;
//
//        if ( data.getColumnDimension() > 1) {
//
//            // compute covariance matrix if we have more than 1 attribute
//            Covariance c = new Covariance(data);
//            covMatrix = c.getCovarianceMatrix();
//
//        } else {
//
//            // if we only have one attribute calculate the variance instead
//            covMatrix = MatrixUtils.createRealMatrix(1,1);
//            covMatrix.setEntry(0, 0, StatUtils.variance(data.getColumn(0)));
//
//        }
//
//        // get the eigenvalues and eigenvectors of the covariance matrixE
//        EigenDecomposition eDecomp = new EigenDecompositionImpl(covMatrix,0.0);
//
//        // set the eigenVectors matrix
//        // the columns of the eigenVectors matrix are the eigenVectors of
//        // the covariance matrix
//        RealMatrix eigenVectors = eDecomp.getV();
//
//        // set the eigenValues vector
////        RealVector eigenValues = new ArrayRealVector(eDecomp.getRealEigenvalues());
//
//        //transform the data
//        RealMatrix pcs = data.multiply(eigenVectors);
//
//        locationPrincipalAxis = pcs.getRow(0);

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
        public final static String VIRUS_LOCATIONS = "virusLocations";
        public final static String SERUM_LOCATIONS = "serumLocations";
        public static final String MDS_DIMENSION = "mdsDimension";
        public static final String MDS_PRECISION = "mdsPrecision";

        public static final String LOG_2_TRANSFORM = "log2Transform";
        public static final String TITRATION_THRESHOLD = "titrationThreshold";

        public String getParserName() {
            return ANTIGENIC_TRAIT_LIKELIHOOD;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            String fileName = xo.getStringAttribute(FILE_NAME);
            DataTable<String[]> assayTable;
            try {
                assayTable = DataTable.Text.parse(new FileReader(fileName));
            } catch (IOException e) {
                throw new XMLParseException("Unable to read assay data from file, " + fileName);
            }

            int mdsDimension = xo.getIntegerAttribute(MDS_DIMENSION);

            boolean log2Transform = false;
            if (xo.hasAttribute(LOG_2_TRANSFORM)) {
                log2Transform = xo.getBooleanAttribute(LOG_2_TRANSFORM);
            }

            double threshold = 0.0;
            if (xo.hasAttribute(TITRATION_THRESHOLD)) {
                threshold = xo.getDoubleAttribute(TITRATION_THRESHOLD);
            }


            // This parameter needs to be linked to the one in the IntegratedMultivariateTreeLikelihood (I suggest that the parameter is created
            // here and then a reference passed to IMTL - which optionally takes the parameter of tip trait values, in which case it listens and
            // updates accordingly.
            CompoundParameter tipTraitParameter = null;
            if (xo.hasChildNamed(TIP_TRAIT)) {
                tipTraitParameter = (CompoundParameter) xo.getElementFirstChild(TIP_TRAIT);
            }

            MatrixParameter virusLocationsParameter = (MatrixParameter) xo.getElementFirstChild(VIRUS_LOCATIONS);
            MatrixParameter serumLocationsParameter = null;

            if (xo.hasChildNamed(SERUM_LOCATIONS)) {
                serumLocationsParameter = (MatrixParameter) xo.getElementFirstChild(SERUM_LOCATIONS);
            }

            if (serumLocationsParameter != null && serumLocationsParameter.getColumnDimension() != virusLocationsParameter.getColumnDimension()) {
                throw new XMLParseException("Virus Locations parameter and Serum Locations parameter have different column dimensions");
            }

            Parameter mdsPrecision = (Parameter) xo.getElementFirstChild(MDS_PRECISION);

            return new AntigenicTraitLikelihood(mdsDimension, mdsPrecision, tipTraitParameter, virusLocationsParameter, serumLocationsParameter, assayTable, log2Transform, threshold);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "Provides the likelihood of immunological assay data such as Hemagglutinin inhibition (HI) given vectors of coordinates" +
                    "for viruses and sera/antisera in some multidimensional 'antigenic' space.";
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newStringRule(FILE_NAME, false, "The name of the file containing the assay table"),
                AttributeRule.newIntegerRule(MDS_DIMENSION, false, "The dimension of the space for MDS"),
                AttributeRule.newBooleanRule(LOG_2_TRANSFORM, true, "Whether to log2 transform the data"),
                AttributeRule.newDoubleRule(TITRATION_THRESHOLD, true, "Titration threshold below which the measurement is not valid"),
                new ElementRule(TIP_TRAIT, CompoundParameter.class, "The parameter of tip locations from the tree", true),
                new ElementRule(VIRUS_LOCATIONS, MatrixParameter.class),
                new ElementRule(SERUM_LOCATIONS, MatrixParameter.class, "An optional set of serum locations", true),
                new ElementRule(MDS_PRECISION, Parameter.class)
        };

        public Class getReturnType() {
            return AntigenicTraitLikelihood.class;
        }
    };

    private final double[][] assayTable;

    private final int tipCount;
    private final int virusLocationCount;
    private final int serumLocationCount;
    private final int[] tipIndices;
    private final int[] virusIndices;

    private final CompoundParameter tipTraitParameter;
    private final MatrixParameter virusLocationsParameter;
    private final MatrixParameter serumLocationsParameter;
    private final Parameter mdsParameter;

    private final int totalMeasurementCount;

    // a set of vectors for each virus giving serum indices for which assay data is available
    private final int[][] measuredSerumIndices;

    private final double titrationThreshold;

    private boolean likelihoodKnown = false;
    private double logLikelihood;
    private double storedLogLikelihood;

    private double[] locationMean;
    private double[] locationPrincipalAxis;

    private boolean distancesKnown = false;
    private double sumOfSquaredResiduals;
    private double storedSumOfSquaredResiduals;
    private double[] distances;
    private double[] storedDistances;

    private boolean[] virusLocationUpdates;
    private boolean[] serumLocationUpdates;
    private boolean[] distanceUpdate;

    private boolean truncationKnown = false;
    private double truncationSum;
    private double storedTruncationSum;
    private double[] truncations;
    private double[] storedTruncations;

    private final boolean isLeftTruncated;
    private final int mdsDimension;
}

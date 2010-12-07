package dr.evomodel.continuous;

import dr.inference.model.*;
import dr.math.MathUtils;
import dr.math.distributions.NormalDistribution;
import dr.util.DataTable;
import dr.xml.*;
import org.apache.commons.math.linear.*;
import org.apache.commons.math.stat.StatUtils;
import org.apache.commons.math.stat.correlation.Covariance;

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
            DataTable<double[]> dataTable,
            String constraint1,
            String constraint2,
            final boolean log2Transform) {

        super(ANTIGENIC_TRAIT_LIKELIHOOD);

        this.mdsDimension = mdsDimension;

        String[] virusNames = dataTable.getRowLabels();
        String[] serumNames = dataTable.getColumnLabels();

//        mdsDimension = virusLocationsParameter.getColumnDimension();

        // the total number of viruses is the number of rows in the table
        virusCount = dataTable.getRowCount();
        // the number of sera is the number of columns
        serumCount = dataTable.getColumnCount();

        tipCount = virusCount;

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
                tipNameMap.put(label, i);

                tipIndices[i] = -1;
            }
        } else {
            tipIndices = null;
        }

        // the virus -> tip map
        virusIndices = new int[virusCount];

        // a set of vectors for each virus giving serum indices for which assay data is available
        measuredSerumIndices = new int[virusCount][];

        // a compressed (no missing values) set of measured assay values between virus and sera.
        this.assayTable = new double[virusCount][];

        int ci1 = -1;
        int ci2 = -1;

        int totalMeasurementCount = 0;
        for (int i = 0; i < virusCount; i++) {
            virusIndices[i] = -1;

            double[] dataRow = dataTable.getRow(i);

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

            if (virusNames[i].equals(constraint1)) {
                ci1 = i;
            }
            if (virusNames[i].equals(constraint2)) {
                ci2 = i;
            }

            int measuredCount = 0;
            for (int j = 0; j < serumCount; j++) {
                if (!Double.isNaN(dataRow[j]) && dataRow[j] > 0) {
                    measuredCount ++;
                }
            }

            assayTable[i] = new double[measuredCount];
            measuredSerumIndices[i] = new int[measuredCount];

            int k = 0;
            for (int j = 0; j < serumCount; j++) {
                if (!Double.isNaN(dataRow[j]) && dataRow[j] > 0) {
                    if (log2Transform) {
                        this.assayTable[i][k] = transform(dataRow[j]);
                    } else {
                        this.assayTable[i][k] = dataRow[j];
                    }
                    measuredSerumIndices[i][k] = j;
                    k ++;
                }
            }
            totalMeasurementCount += measuredCount;
        }

        constraintIndex1 = ci1;
        constraintIndex2 = ci2;

        this.totalMeasurementCount = totalMeasurementCount;

        // a cache of virus to serum distances (serum indices given by array above).
        distances = new double[totalMeasurementCount];
        storedDistances = new double[totalMeasurementCount];

        virusUpdates = new boolean[virusCount];
        serumUpdates = new boolean[serumCount];
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
        virusLocationsParameter.setRowDimension(virusCount);
        addVariable(virusLocationsParameter);

        // some random initial locations
        for (int i = 0; i < virusCount; i++) {
            virusLocationsParameter.getParameter(i).setId(virusNames[i]);
            for (int j = 0; j < mdsDimension; j++) {
                double r = MathUtils.nextGaussian();
                virusLocationsParameter.getParameter(i).setParameterValue(j, r);
            }
        }

        this.serumLocationsParameter = serumLocationsParameter;
        if (virusLocationsParameter != serumLocationsParameter) {
            serumLocationsParameter.setColumnDimension(mdsDimension);
            serumLocationsParameter.setRowDimension(serumCount);
            addVariable(serumLocationsParameter);
        }

        // some random initial locations
        for (int i = 0; i < serumCount; i++) {
            virusLocationsParameter.getParameter(i).setId(serumNames[i]);
            for (int j = 0; j < mdsDimension; j++) {
                double r = MathUtils.nextGaussian();
                serumLocationsParameter.getParameter(i).setParameterValue(j, r);
            }
        }

        this.mdsParameter = mdsPrecision;
        addVariable(mdsPrecision);
   
        this.isLeftTruncated = true; // Re-normalize likelihood for strictly positive distances

        addStatistic(pcaStatistic);
    }

    private double transform(final double value) {
        // transform to log_2
        return Math.log(value) / Math.log(2.0);
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Variable.ChangeType type) {
        // TODO Flag which cachedDistances or mdsPrecision need to be updated

        if (variable == virusLocationsParameter) {
            if (tipTraitParameter != null) {
                // the virus locations have changed so update the tipTraitParameter
                int k = 0;
                for (int i = 0; i < tipCount; i++) {
                    if (tipIndices[i] != -1) {
                        Parameter virusLoc = virusLocationsParameter.getParameter(tipIndices[i]);
                        for (int j = 0; j < mdsDimension; j++) {
                            tipTraitParameter.setParameterValue(k, virusLoc.getValue(j));
                            k++;
                        }
                    } else {
                        k += mdsDimension;
                    }
                }
            }

            if (constraintIndex1 >= 0) {
                virusLocationsParameter.getParameter(constraintIndex1).setParameterValueQuietly(0, 0.0);
                virusLocationsParameter.getParameter(constraintIndex1).setParameterValueQuietly(1, 0.0);
            }
            if (constraintIndex2 >= 0) {
                virusLocationsParameter.getParameter(constraintIndex2).setParameterValueQuietly(0, 0.0);
            }

            virusUpdates[index / mdsDimension] = true;
            distancesKnown = false;
        } else if (variable == serumLocationsParameter) {
            serumUpdates[index / mdsDimension] = true;
            distancesKnown = false;
        } else if (variable == mdsParameter) {
            for (int i = 0; i < distanceUpdate.length; i++) {
                distanceUpdate[i] = true;
            }
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

        for (int i = 0; i < virusUpdates.length; i++) {
            virusUpdates[i] = true;
        }

        for (int i = 0; i < serumUpdates.length; i++) {
            serumUpdates[i] = true;
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

        for (int i = 0; i < virusUpdates.length; i++) {
            virusUpdates[i] = false;
        }

        for (int i = 0; i < serumUpdates.length; i++) {
            serumUpdates[i] = false;
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
                if (virusUpdates[i] || serumUpdates[measuredSerumIndices[i][j]]) {
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

    private double[][] calculatePCA() {
        double[][] locations = new double[virusCount + serumCount][mdsDimension];
        int k = 0;
        for (int i = 0; i < virusCount; i++) {
            for (int j = 0; j < mdsDimension; j++) {
                locations[k][j] = virusLocationsParameter.getParameter(i).getParameterValue(j);
            }
            k++;
        }
        for (int i = 0; i < serumCount; i++) {
            for (int j = 0; j < mdsDimension; j++) {
                locations[k][j] = serumLocationsParameter.getParameter(i).getParameterValue(j);
            }
            k++;
        }

        RealMatrix data = MatrixUtils.createRealMatrix(locations);
        // compute the covariance matrix
        RealMatrix covMatrix = null;

        if ( data.getColumnDimension() > 1) {

            // compute covariance matrix if we have more than 1 attribute
            Covariance c = new Covariance(data);
            covMatrix = c.getCovarianceMatrix();

        } else {

            // if we only have one attribute calculate the variance instead
            covMatrix = MatrixUtils.createRealMatrix(1,1);
            covMatrix.setEntry(0, 0, StatUtils.variance(data.getColumn(0)));

        }

        // get the eigenvalues and eigenvectors of the covariance matrixE
        EigenDecomposition eDecomp = new EigenDecompositionImpl(covMatrix,0.0);

        // set the eigenVectors matrix
        // the columns of the eigenVectors matrix are the eigenVectors of
        // the covariance matrix
        RealMatrix eigenVectors = eDecomp.getV();

        // set the eigenValues vector
//        RealVector eigenValues = new ArrayRealVector(eDecomp.getRealEigenvalues());

        //transform the data
        RealMatrix pcs = data.multiply(eigenVectors);

        return pcs.getData();
    }

    private final Statistic pcaStatistic = new Statistic.Abstract() {

        public String getStatisticName() {
            return "pca";
        }

        public int getDimension() {
            return 1;
        }

        public double getStatisticValue(int dim) {
            double[][] pcs  = calculatePCA();
            return pcs[0][0];
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
        public static final String CONSTRAINT_1 = "constraint1";
        public static final String CONSTRAINT_2 = "constraint2";

        public static final String LOG_2_TRANSFORM = "log2Transform";

        public String getParserName() {
            return ANTIGENIC_TRAIT_LIKELIHOOD;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            String fileName = xo.getStringAttribute(FILE_NAME);
            DataTable<double[]> assayTable;
            try {
                assayTable = DataTable.Double.parse(new FileReader(fileName));
            } catch (IOException e) {
                throw new XMLParseException("Unable to read assay data from file, " + fileName);
            }

            int mdsDimension = xo.getIntegerAttribute(MDS_DIMENSION);

            boolean log2Transform = false;
            if (xo.hasAttribute(LOG_2_TRANSFORM)) {
                xo.getBooleanAttribute(LOG_2_TRANSFORM);
            }

            String constraint1 = null;
            if (xo.hasAttribute(CONSTRAINT_1)) {
                constraint1 = xo.getStringAttribute(CONSTRAINT_1);
            }
            String constraint2 = null;
            if (xo.hasAttribute(CONSTRAINT_2)) {
                constraint2 = xo.getStringAttribute(CONSTRAINT_2);
            }

            // This parameter needs to be linked to the one in the IntegratedMultivariateTreeLikelihood (I suggest that the parameter is created
            // here and then a reference passed to IMTL - which optionally takes the parameter of tip trait values, in which case it listens and
            // updates accordingly.
            CompoundParameter tipTraitParameter = null;
            if (xo.hasChildNamed(TIP_TRAIT)) {
                tipTraitParameter = (CompoundParameter) xo.getElementFirstChild(TIP_TRAIT);
            }

            MatrixParameter virusLocationsParameter = (MatrixParameter) xo.getElementFirstChild(VIRUS_LOCATIONS);
            MatrixParameter serumLocationsParameter = (MatrixParameter) xo.getElementFirstChild(SERUM_LOCATIONS);

            if (serumLocationsParameter.getColumnDimension() != virusLocationsParameter.getColumnDimension()) {
                throw new XMLParseException("Virus Locations parameter and Serum Locations parameter have different column dimensions");
            }

            Parameter mdsPrecision = (Parameter) xo.getElementFirstChild(MDS_PRECISION);

            return new AntigenicTraitLikelihood(mdsDimension, mdsPrecision, tipTraitParameter, virusLocationsParameter, serumLocationsParameter, assayTable, constraint1, constraint2, log2Transform);
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
                AttributeRule.newStringRule(CONSTRAINT_1, true, "The name of a virus to act as a constraint"),
                AttributeRule.newStringRule(CONSTRAINT_2, true, "The name of a virus to act as a constraint"),
                new ElementRule(TIP_TRAIT, CompoundParameter.class, "The parameter of tip locations from the tree", true),
                new ElementRule(VIRUS_LOCATIONS, MatrixParameter.class),
                new ElementRule(SERUM_LOCATIONS, MatrixParameter.class),
                new ElementRule(MDS_PRECISION, Parameter.class)
        };

        public Class getReturnType() {
            return AntigenicTraitLikelihood.class;
        }
    };

    public static void main(String[] argv) {
        DataTable<double[]> assayTable = null;
        try {
            assayTable = DataTable.Double.parse(new FileReader(argv[0]));
        } catch (IOException e) {
            System.err.println("Unable to read assay data from file, " + argv[0]);
            System.exit(1);
        }

        String[] virusNames = assayTable.getRowLabels();
        String[] serumNames = assayTable.getColumnLabels();

        Set<String> serumSet = new HashSet<String>();
        Map<String, String> serumMap = new HashMap<String, String>();
        Map<String, String> virusMap = new HashMap<String, String>();

        for (String serum : serumNames) {
            String[] bits = serum.split("/");
            String country = bits[0];

            String name = "";
            if (bits.length > 1) {
                bits = bits[1].split("-");
                String code = bits[0];
                name = country + "_" + code;
            } else {
                if (country.startsWith("IVR-134")) {
                    name = "WYOMING_3";
                } else if (country.startsWith("IVR-135")) {
                    name = "KUMAMOTO_102";
                } else if (country.startsWith("RESVIR")) {
                    name = "PANAMA_2007";
                } else {
                    throw new RuntimeException("Unknown virus");
                }
            }

            serumSet.add(name);
            serumMap.put(serum, name);

            for (String virus : virusNames) {
                if (virus.toUpperCase().contains(name + "_")) {
                    virusMap.put(virus, name);
                }
            }
        }

        Map<String, double[]> columnMap = new HashMap<String, double[]>();

        for (int j = 0; j < serumNames.length; j++) {
            String serum = serumNames[j];

            double[] column = assayTable.getColumn(j);
            String code = serumMap.get(serum);

            double[] column0 = columnMap.get(code);
            if (column0 != null) {
                for (int i = 0; i < virusNames.length; i++) {
                    if (!Double.isNaN(column[i]) && !Double.isNaN(column0[i])) {
                        System.out.println("Duplicate - " + serum + ": " + column[i] +  ", " + column0[i]);
                    }
                }
            } else {
                columnMap.put(code, column);
            }


        }

        System.out.println("Unique serum count: " + serumSet.size());
        System.out.println("Unique serum+virus count: " + virusMap.keySet().size());
    }

    private final double[][] assayTable;

    private final int tipCount;
    private final int virusCount;
    private final int serumCount;
    private final int[] tipIndices;
    private final int[] virusIndices;

    private final int constraintIndex1;
    private final int constraintIndex2;

    private final CompoundParameter tipTraitParameter;
    private final MatrixParameter virusLocationsParameter;
    private final MatrixParameter serumLocationsParameter;
    private final Parameter mdsParameter;

    private final int totalMeasurementCount;

    // a set of vectors for each virus giving serum indices for which assay data is available
    private final int[][] measuredSerumIndices;


    private boolean likelihoodKnown = false;
    private double logLikelihood;
    private double storedLogLikelihood;

    private boolean distancesKnown = false;
    private double sumOfSquaredResiduals;
    private double storedSumOfSquaredResiduals;
    private double[] distances;
    private double[] storedDistances;

    private boolean[] virusUpdates;
    private boolean[] serumUpdates;
    private boolean[] distanceUpdate;

    private boolean truncationKnown = false;
    private double truncationSum;
    private double storedTruncationSum;
    private double[] truncations;
    private double[] storedTruncations;

    private final boolean isLeftTruncated;
    private final int mdsDimension;
}

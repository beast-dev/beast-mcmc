package dr.evomodel.continuous;

import dr.inference.model.*;
import dr.math.distributions.NormalDistribution;
import dr.util.DataTable;
import dr.xml.*;

import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

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
            DataTable<double[]> dataTable) {

        super(ANTIGENIC_TRAIT_LIKELIHOOD);

        this.mdsDimension = mdsDimension;

        String[] virusNames = dataTable.getRowLabels();
        String[] serumNames = dataTable.getColumnLabels();

//        mdsDimension = virusLocationsParameter.getColumnDimension();

        // the total number of viruses is the number of rows in the table
        int virusCount = dataTable.getRowCount();
        // the number of sera is the number of columns
        int serumCount = dataTable.getColumnCount();

        tipCount = virusCount;

        Map<String, Integer> tipNameMap = null;
        if (tipTraitParameter != null) {
            if (tipCount != tipTraitParameter.getNumberOfParameters()) {
                System.err.println("Tree has different number of tips than the number of viruses");
            }

            // the tip -> virus map
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
                    this.assayTable[i][k] = transform(dataRow[j]); // TODO Code review here, was dataRow[k]
                    measuredSerumIndices[i][k] = j;
                    k ++;
                }
            }
            totalMeasurementCount += measuredCount;
        }

        this.totalMeasurementCount = totalMeasurementCount;

        // a cache of virus to serum distances (serum indices given by array above).
        distances = new double[totalMeasurementCount];
        storedDistances = new double[totalMeasurementCount];

        // a cache of individual truncations
//        truncations = new double[totalMeasurementCount];
//        storedTruncations = new double[totalMeasurementCount];

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

        this.serumLocationsParameter = serumLocationsParameter;
        serumLocationsParameter.setColumnDimension(mdsDimension);
        serumLocationsParameter.setRowDimension(serumCount);
        addVariable(serumLocationsParameter);


        this.mdsParameter = mdsPrecision;
        addVariable(mdsPrecision);

        this.isLeftTruncated = true; // Re-normalize likelihood for strictly positive distances
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
            
            distancesKnown = false;
        } else  if (variable == serumLocationsParameter) {
            distancesKnown = false;
        }
        distancesKnown = false;
        truncationKnown = false;
        likelihoodKnown = false;
    }

    @Override
    protected void storeState() {
//        System.arraycopy(distances, 0, storedDistances, 0, distances.length);
//        System.arraycopy(truncations, 0, storedTruncations, 0, truncations.length);

        storedLogLikelihood = logLikelihood;
        storedTruncation = truncation;
        storedSumOfSquaredResiduals = sumOfSquaredResiduals;
    }

    @Override
    protected void restoreState() {
//        double[] tmp = storedDistances;
//        storedDistances = distances;
//        distances = tmp;
//        distancesKnown = true;

//        tmp = storedTruncations;
//        storedTruncations = truncations;
//        truncations = tmp;

        logLikelihood = storedLogLikelihood;
        likelihoodKnown = true;

        truncation = storedTruncation;
        truncationKnown = true;

        sumOfSquaredResiduals = storedSumOfSquaredResiduals;

        makeDirty();
    }

    @Override
    protected void acceptState() {
        // do nothing
    }

    public void makeDirty() {
        distancesKnown = false;
        likelihoodKnown = false;
        truncationKnown = false;
    }

    public Model getModel() {
        return this;
    }

    public double getLogLikelihood() {
        // TODO Only recompute when not known: distances or mdsPrecision changed
        if (!likelihoodKnown) {
            if (!distancesKnown) {
                calculateDistances();
                sumOfSquaredResiduals = calculateSumOfSquaredResiduals();
                distancesKnown = true;
            }

            logLikelihood = computeLogLikelihood();
            likelihoodKnown = true;
        }

        return logLikelihood;
    }

    // This function can be overwritten to implement other sampling densities, i.e. discrete ranks
    protected double computeLogLikelihood() {

        double precision = mdsParameter.getParameterValue(0);
        double logLikelihood = (totalMeasurementCount / 2) * Math.log(precision) - 0.5 * precision * sumOfSquaredResiduals;

        if (isLeftTruncated) {
            if (!truncationKnown) {
                truncation = calculateTruncation(precision);
                truncationKnown = true;
            }
            logLikelihood -= truncation;
        }
        return logLikelihood;
    }

    private double calculateTruncation(double precision) {
        double sum = 0.0;
        double sd = 1.0 / Math.sqrt(precision);
        int k = 0;
        for (int i = 0; i < assayTable.length; i++) {
            for (int j = 0; j < assayTable[i].length; j++) {
                double t = Math.log(NormalDistribution.cdf(distances[k], 0.0, sd));
//                truncations[k] = t;
                sum += t;
                k++;
            }
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
                distances[k] = calculateDistance(virusLocationsParameter.getParameter(i),
                        serumLocationsParameter.getParameter(measuredSerumIndices[i][j]));
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

            return new AntigenicTraitLikelihood(mdsDimension, mdsPrecision, tipTraitParameter, virusLocationsParameter, serumLocationsParameter, assayTable);
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
                new ElementRule(TIP_TRAIT, CompoundParameter.class, "The parameter of tip locations from the tree", true),
                new ElementRule(VIRUS_LOCATIONS, MatrixParameter.class),
                new ElementRule(SERUM_LOCATIONS, MatrixParameter.class),
                new ElementRule(MDS_PRECISION, Parameter.class)
        };

        public Class getReturnType() {
            return AntigenicTraitLikelihood.class;
        }
    };

    private final double[][] assayTable;

    private final int tipCount;
    private final int[] tipIndices;
    private final int[] virusIndices;

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

    private boolean truncationKnown = false;
    private double truncation;
    private double storedTruncation;
//    private double[] truncations;
//    private double[] storedTruncations;


    private final boolean isLeftTruncated;
    private final int mdsDimension;
}

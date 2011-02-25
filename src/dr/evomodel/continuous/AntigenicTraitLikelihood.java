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
public class AntigenicTraitLikelihood extends MultidimensionalScalingLikelihood  {

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

        String[] virusNames = dataTable.getRowLabels();
        boolean hasAliases = false;
        if (virusNames[0].equals("aliases")) {
            // if the first row is labelled aliases then there are aliases so the
            // viruses start on row 2.
            hasAliases = true;

            // Java 1.6:
//            virusNames = Arrays.copyOfRange(virusNames, 1, virusNames.length);

            // Java 1.5:
            virusNames = new String[virusNames.length - 1];
            for (int i = 0; i < virusNames.length; i++) {
                virusNames[i] = dataTable.getRowLabels()[i + 1];
            }

        }

        // the total number of viruses is the number of rows in the table
        int virusCount = dataTable.getRowCount() - (hasAliases ? 1 : 0);
        int assayCount = dataTable.getColumnCount();
        int serumCount;

        String[] assayNames = dataTable.getColumnLabels();
        int[] assayToSerumIndices = new int[assayNames.length];

        Set<String> aliasSet = null;
        List<String> aliasNames = null;

        String[] serumNames = null;

        if (hasAliases) {
            aliasSet = new HashSet<String>();
            aliasNames = new ArrayList<String>();

            String[] aliases = dataTable.getRow(0);
            for (int i = 0; i < assayNames.length; i++) {
                if (aliasSet.contains(aliases[i])) {
                    assayToSerumIndices[i] = aliasNames.indexOf(aliases[i]);
                } else {
                    aliasSet.add(aliases[i]);
                    aliasNames.add(aliases[i]);
                    assayToSerumIndices[i] = aliasNames.size() - 1;
                }
            }

            // the number of serum locations is the number of aliases
            serumCount = aliasNames.size();
            serumNames = new String[aliasNames.size()];
            aliasNames.toArray(serumNames);

        } else {
            // the number of serum locations is the number of columns
            serumCount = assayCount;

            // one alias for one serum
            for (int i = 0; i < assayToSerumIndices.length; i++) {
                assayToSerumIndices[i] = i;
            }

            serumNames = assayNames;
        }


        int start = hasAliases ? 1 : 0;

        List<Double> observationList = new ArrayList<Double>();
        List<Integer> distanceIndexList = new ArrayList<Integer>();
        List<Integer> rowIndexList = new ArrayList<Integer>();
        List<Integer> columnIndexList = new ArrayList<Integer>();
        List<ObservationType> observationTypeList = new ArrayList<ObservationType>();

        int[] virusObservationCounts = new int[virusCount];
        int[] serumObservationCounts = new int[serumCount];

        // the largest measured value for any given column of data
        // Currently this is the largest across any assay column for a given antisera.
        // Optionally could normalize by individual assay column
        double[] maxAssayValue = new double[serumCount];

        // Build a sparse matrix of non-missing assay values
        int u = 0;
        for (int i = 0; i < virusCount; i++) {
            String[] dataRow = dataTable.getRow(i + start);

            for (int j = 0; j < serumCount; j++) {

                boolean first = true;

                for (int k = 0; k < assayCount; k++) {
                    if (assayToSerumIndices[k] == j) {
                        Double value = null;
                        ObservationType type = null;

                        if (dataRow[k].startsWith("<")) {
                            // is a lower bound
                            value = convertString(dataRow[k].substring(1));
                            if (Double.isNaN(value)) {
                                throw new RuntimeException("Illegal value in table as a threshold");
                            }
                            type = ObservationType.LOWER_BOUND;
                        } else if (dataRow[k].startsWith(">")) {
                            // is a lower bound
                            value = convertString(dataRow[k].substring(1));
                            if (Double.isNaN(value)) {
                                throw new RuntimeException("Illegal value in table as a threshold");
                            }
                            type = ObservationType.UPPER_BOUND;
                        } else {
                            value = convertString(dataRow[k]);
                            type = ObservationType.POINT;
                        }
                        if (!Double.isNaN(value)) {
                            observationList.add(value);
                            observationTypeList.add(type);
                            distanceIndexList.add(u);
                            virusObservationCounts[i]++;
                            serumObservationCounts[j]++;

                            if (value > maxAssayValue[j]) {
                                maxAssayValue[j] = value;
                            }
                        }

                        if (first) {
                            // if this is the first time an observation for this virus/serum pair is found:
                            rowIndexList.add(i);
                            columnIndexList.add(j);
                            first = false;
                            u++;
                        }

                    }
                }

            }
        }

        // check that all the viruses and sera have observations
        for (int i = 0; i < virusCount; i++) {
            if (virusObservationCounts[i] == 0) {
                System.err.println("WARNING: Virus " + virusNames[i] + " has 0 observations");
            }
        }
        for (int j = 0; j < serumCount; j++) {
            if (serumObservationCounts[j] == 0) {
                System.err.println("WARNING: Antisera " + serumNames[j] + " has 0 observations");
            }
        }


        // Convert into arrays
        double[] observations = new double[observationList.size()];
        for (int i = 0; i < observationList.size(); i++) {
            observations[i] = observationList.get(i);
        }

        int[] distanceIndices = new int[distanceIndexList.size()];
        for (int i = 0; i < distanceIndexList.size(); i++) {
            distanceIndices[i] = distanceIndexList.get(i);
        }

        int[] rowIndices = new int[rowIndexList.size()];
        for (int i = 0; i < rowIndexList.size(); i++) {
            rowIndices[i] = rowIndexList.get(i);
        }

        int[] columnIndices = new int[columnIndexList.size()];
        for (int i = 0; i < columnIndexList.size(); i++) {
            columnIndices[i] = columnIndexList.get(i);
        }

        ObservationType[] observationTypes = new ObservationType[observationTypeList.size()];
        observationTypeList.toArray(observationTypes);

        // transform and normalize the data if required
        if (log2Transform) {
            for (int i = 0; i < observations.length; i++) {
                observations[i] = transform(observations[i], maxAssayValue[columnIndices[i]]);
                // the transformation reverses the bounds
                observationTypes[i] = (observationTypes[i] == ObservationType.UPPER_BOUND ? ObservationType.LOWER_BOUND : (observationTypes[i] == ObservationType.LOWER_BOUND ? ObservationType.UPPER_BOUND : observationTypes[i]));
            }
        }

        initialize(mdsDimension, mdsPrecision, tipTraitParameter, virusLocationsParameter, serumLocationsParameter, virusNames, serumNames, observations, observationTypes, distanceIndices, rowIndices, columnIndices);

        // some random initial locations
//        for (int i = 0; i < virusCount; i++) {
//            virusLocationsParameter.getParameter(i).setId(virusNames[i]);
//            for (int j = 0; j < mdsDimension; j++) {
//                double r = MathUtils.nextGaussian();
//                virusLocationsParameter.getParameter(i).setParameterValue(j, r);
//            }
//        }

        // some random initial locations
//            for (int i = 0; i < serumCount; i++) {
//                serumLocationsParameter.getParameter(i).setId(serumNames[i]);
//                for (int j = 0; j < mdsDimension; j++) {
//                    double r = MathUtils.nextGaussian();
//                    serumLocationsParameter.getParameter(i).setParameterValue(j, r);
//                }
//            }
    }

    private double convertString(String value) {
        try {
            return java.lang.Double.valueOf(value);
        } catch (NumberFormatException nfe) {
            return java.lang.Double.NaN;
        }
    }

    private double transform(final double value, final double maxValue) {
        return Math.log(maxValue / value) / Math.log(2.0);
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

}

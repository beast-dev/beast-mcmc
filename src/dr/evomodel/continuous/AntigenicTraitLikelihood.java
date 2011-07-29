package dr.evomodel.continuous;

import dr.inference.model.*;
import dr.math.MathUtils;
import dr.math.distributions.NormalDistribution;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;
import dr.util.DataTable;
import dr.xml.*;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;
import java.util.logging.Logger;

/**
 * @author Andrew Rambaut
 * @author Marc Suchard
 * @version $Id$
 */
public class AntigenicTraitLikelihood extends MultidimensionalScalingLikelihood implements Citable {

    public final static String ANTIGENIC_TRAIT_LIKELIHOOD = "antigenicTraitLikelihood";

    /**
     * Constructor
     * @param mdsDimension dimension of the mds space
     * @param mdsPrecision parameter which gives the precision of the bmds
     * @param tipTraitParameter a parameter of locations for the tips of the tree (mapped to virus locations)
     * @param locationsParameter a parameter of locations of viruses/sera
     * @param dataTable the assay table (virus in rows, serum assays in columns)
     * @param virusAntiserumMap a map of viruses to corresponding sera
     * @param assayAntiserumMap a map of repeated assays for a given sera
     * @param log2Transform transform the data into log 2 space
     */
    public AntigenicTraitLikelihood(
            int mdsDimension,
            Parameter mdsPrecision,
            CompoundParameter tipTraitParameter,
            MatrixParameter locationsParameter,
            DataTable<String[]> dataTable,
            Map<String, String> virusAntiserumMap,
            Map<String, String> assayAntiserumMap,
            final boolean log2Transform) {

        super(ANTIGENIC_TRAIT_LIKELIHOOD);

        Logger.getLogger("dr.evomodel").info("Using EvolutionaryCartography model. Please cite:\n" + Citable.Utils.getCitationString(this));

        String[] virusNames = dataTable.getRowLabels();

        // the total number of viruses is the number of rows in the table
        int virusCount = dataTable.getRowCount();
        int assayCount = dataTable.getColumnCount();
        int serumCount;

        String[] columnNames = dataTable.getColumnLabels();
        int[] columnToSerumIndices = new int[columnNames.length];

        double[][] observationValueTable = new double[virusCount][assayCount];
        ObservationType[][] observationTypeTable = new ObservationType[virusCount][assayCount];

        initalizeTable(dataTable, observationValueTable, observationTypeTable, log2Transform);

        List<String> aliasNames = null;

        String[] serumNames = null;

        if (assayAntiserumMap != null) {
            aliasNames = new ArrayList<String>(new TreeSet<String>(assayAntiserumMap.values()));

            for (int i = 0; i < columnNames.length; i++) {
                String alias = assayAntiserumMap.get(columnNames[i]);
                if (alias == null) {
                    throw new IllegalArgumentException("Alias missing for assay, " + columnNames[i]);
                }

                columnToSerumIndices[i] = aliasNames.indexOf(alias);
            }

            // the number of serum locations is the number of aliases
            serumCount = aliasNames.size();
            serumNames = new String[aliasNames.size()];
            aliasNames.toArray(serumNames);

        } else {
            // the number of serum locations is the number of columns
            serumCount = assayCount;

            // one alias for one serum
            for (int i = 0; i < columnToSerumIndices.length; i++) {
                columnToSerumIndices[i] = i;
            }

            serumNames = columnNames;
        }

        int[] virusToSerumIndices = new int[virusNames.length];

        if (virusAntiserumMap != null) {
            aliasNames = new ArrayList<String>(new TreeSet<String>(assayAntiserumMap.values()));

            for (int i = 0; i < virusNames.length; i++) {
                String alias = virusAntiserumMap.get(virusNames[i]);
                if (alias != null) {
                    virusToSerumIndices[i] = aliasNames.indexOf(alias);
                } else {
                    virusToSerumIndices[i] = -1;
                }
            }

        } else {

            // try and match rows to columns using their names
            List<String> serumNameList = Arrays.asList(serumNames);
            for (int i = 0; i < virusToSerumIndices.length; i++) {
                virusToSerumIndices[i] = serumNameList.indexOf(virusNames[i]);
                if (virusToSerumIndices[i] == -1) {
                    throw new IllegalArgumentException("Row, " + virusNames[i] + ", doesn't match any column heading");
                }
            }
        }

        List<Double> observationList = new ArrayList<Double>();
        List<Integer> distanceIndexList = new ArrayList<Integer>();
        List<Integer> rowIndexList = new ArrayList<Integer>();
        List<Integer> columnIndexList = new ArrayList<Integer>();
        List<ObservationType> observationTypeList = new ArrayList<ObservationType>();

        int[] virusObservationCounts = new int[virusCount];
        int[] serumObservationCounts = new int[serumCount];

        int virusSerumPairObservationCounts = 0;


        // Build a sparse matrix of non-missing assay values
        int u = 0;
        for (int i = 0; i < virusCount; i++) {

            for (int j = 0; j < serumCount; j++) {

                boolean isVirusSerumPair = (virusToSerumIndices[i] == j);

                boolean first = true;

                for (int k = 0; k < assayCount; k++) {
                    if (columnToSerumIndices[k] == j) {
                        Double value = observationValueTable[i][k];
                        ObservationType type = observationTypeTable[i][k];


                        if (type != ObservationType.MISSING) {
                            observationList.add(value);
                            observationTypeList.add(type);
                            distanceIndexList.add(u);
                            virusObservationCounts[i]++;
                            serumObservationCounts[j]++;

                            if(isVirusSerumPair) {
                                virusSerumPairObservationCounts ++;
                            }

                            rowIndexList.add(i);
                            columnIndexList.add(j);
                        }

                        if (first) {
                            // if this is the first time an observation for this virus/serum pair is found:
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

        // **** not finished - need to handle not matched virus/serum
        int[] rowLocationIndices = new int[rowIndexList.size()];
        for (int i = 0; i < rowLocationIndices.length; i++) {
            rowLocationIndices[i] = i;
        }
        int[] columnLocationIndices = new int[columnIndexList.size()];
        for (int i = 0; i < columnLocationIndices.length; i++) {
            columnLocationIndices[i] = i;
        }

        ObservationType[] observationTypes = new ObservationType[observationTypeList.size()];
        observationTypeList.toArray(observationTypes);

        int thresholdCount = 0;
        for (int i = 0; i < observations.length; i++) {
            thresholdCount += (observationTypes[i] != ObservationType.POINT ? 1 : 0);
        }

        StringBuilder sb = new StringBuilder();
        sb.append("\tAntigenicTraitLikelihood:\n");
        sb.append("\t\t" + virusNames.length + " viruses\n");
        sb.append("\t\t" + serumNames.length + " antisera\n");
        sb.append("\t\t" + observations.length + " observations\n");
        sb.append("\t\t" + virusSerumPairObservationCounts + " observations for virus/serum pairs\n");
        sb.append("\t\t" + thresholdCount + " threshold observations\n");
        Logger.getLogger("dr.evomodel").info(sb.toString());

        initialize(mdsDimension, mdsPrecision, tipTraitParameter, locationsParameter, virusNames, serumNames, observations, observationTypes, distanceIndices, rowLocationIndices, columnLocationIndices);

        // some random initial locations
        for (int i = 0; i < virusCount; i++) {
            locationsParameter.getParameter(i).setId(virusNames[i]);
            for (int j = 0; j < mdsDimension; j++) {
                double r = MathUtils.nextGaussian();
                locationsParameter.getParameter(i).setParameterValue(j, r);

//                locationsParameter.getParameter(i).setParameterValue(j, i * 1000);
            }
        }


    }

    private void initalizeTable(DataTable<String[]> dataTable, double[][] observationValueTable, ObservationType[][] observationTypeTable, boolean log2Transform) {
        // the largest measured value for any given column of data
        double[] maxColumnValue = new double[dataTable.getColumnCount()];
        // the largest measured value over all
        double maxAssayValue = 0.0;

        for (int i = 0; i < dataTable.getRowCount(); i++) {
            String[] dataRow = dataTable.getRow(i);

            for (int j = 0; j < dataTable.getColumnCount(); j++) {
                Double value = null;
                ObservationType type = null;

                if (dataRow[j].startsWith("<")) {
                    // is a lower bound
                    value = convertString(dataRow[j].substring(1));
                    if (Double.isNaN(value)) {
                        throw new RuntimeException("Illegal value in table as a threshold");
                    }
                    type = ObservationType.LOWER_BOUND;
                } else if (dataRow[j].startsWith(">")) {
                    // is a lower bound
                    value = convertString(dataRow[j].substring(1));
                    if (Double.isNaN(value)) {
                        throw new RuntimeException("Illegal value in table as a threshold");
                    }
                    type = ObservationType.UPPER_BOUND;
                } else {
                    value = convertString(dataRow[j]);
                    type = Double.isNaN(value) ? ObservationType.MISSING : ObservationType.POINT;
                }

                observationValueTable[i][j] = value;
                observationTypeTable[i][j] = type;

                if (!Double.isNaN(value)) {
                    if (value > maxColumnValue[j]) {
                        maxColumnValue[j] = value;
                    }
                    if (value > maxAssayValue) {
                        maxAssayValue = value;
                    }

                }

            }
        }

        if (log2Transform) {
            // transform and normalize the data...
            for (int i = 0; i < dataTable.getRowCount(); i++) {
                for (int j = 0; j < dataTable.getColumnCount(); j++) {
                    observationValueTable[i][j] = transform(observationValueTable[i][j], maxColumnValue[j], 2);
                    // the transformation reverses the bounds
                    if (observationTypeTable[i][j] == ObservationType.UPPER_BOUND) {
                        observationTypeTable[i][j] = ObservationType.LOWER_BOUND;
                    } else if (observationTypeTable[i][j] == ObservationType.LOWER_BOUND) {
                        observationTypeTable[i][j] = ObservationType.UPPER_BOUND;
                    }
                }
            }
        }
    }

    private double convertString(String value) {
        try {
            return java.lang.Double.valueOf(value);
        } catch (NumberFormatException nfe) {
            return java.lang.Double.NaN;
        }
    }

    private double transform(final double value, final double maxValue, final double base) {
        // log2(maxValue / value)
        return (Math.log(maxValue) - Math.log(value)) / Math.log(base);
    }

    // **************************************************************
    // XMLObjectParser
    // **************************************************************

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {
        public final static String FILE_NAME = "fileName";
        public final static String VIRUS_MAP_FILE_NAME = "virusMapFile";
        public final static String ASSAY_MAP_FILE_NAME = "assayMapFile";

        public final static String TIP_TRAIT = "tipTrait";
        public final static String LOCATIONS = "locations";
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

            Map<String, String> virusAntiserumMap = null;
            if (xo.hasAttribute(VIRUS_MAP_FILE_NAME)) {
                try {
                    virusAntiserumMap = readMap(xo.getStringAttribute(VIRUS_MAP_FILE_NAME));
                } catch (IOException e) {
                    throw new XMLParseException("Virus map file not found: " + xo.getStringAttribute(VIRUS_MAP_FILE_NAME));
                }
            }

            Map<String, String> assayAntiserumMap = null;
            if (xo.hasAttribute(ASSAY_MAP_FILE_NAME)) {
                try {
                    assayAntiserumMap = readMap(xo.getStringAttribute(ASSAY_MAP_FILE_NAME));
                } catch (IOException e) {
                    throw new XMLParseException("Assay map file not found: " + xo.getStringAttribute(ASSAY_MAP_FILE_NAME));
                }
            }

            int mdsDimension = xo.getIntegerAttribute(MDS_DIMENSION);

            boolean log2Transform = false;
            if (xo.hasAttribute(LOG_2_TRANSFORM)) {
                log2Transform = xo.getBooleanAttribute(LOG_2_TRANSFORM);
            }

            // This parameter needs to be linked to the one in the IntegratedMultivariateTreeLikelihood (I suggest that the parameter is created
            // here and then a reference passed to IMTL - which optionally takes the parameter of tip trait values, in which case it listens and
            // updates accordingly.
            CompoundParameter tipTraitParameter = null;
            if (xo.hasChildNamed(TIP_TRAIT)) {
                tipTraitParameter = (CompoundParameter) xo.getElementFirstChild(TIP_TRAIT);
            }

            MatrixParameter locationsParameter = (MatrixParameter) xo.getElementFirstChild(LOCATIONS);

            Parameter mdsPrecision = (Parameter) xo.getElementFirstChild(MDS_PRECISION);

            return new AntigenicTraitLikelihood(mdsDimension, mdsPrecision, tipTraitParameter, locationsParameter, assayTable, virusAntiserumMap, assayAntiserumMap, log2Transform);
        }

        private  Map<String, String> readMap(String fileName) throws IOException {
            BufferedReader reader = new BufferedReader(new FileReader(fileName));

            Map<String, String> map = new HashMap<String, String>();

            String line = reader.readLine();
            while (line != null) {
                if (line.trim().length() > 0) {
                    String[] parts = line.split("\t");
                    if (parts.length > 1) {
                        map.put(parts[0], parts[1]);
                    }
                }
                line = reader.readLine();
            }

            reader.close();

            return map;
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
                AttributeRule.newStringRule(VIRUS_MAP_FILE_NAME, true, "The name of the file containing the virus to serum map"),
                AttributeRule.newStringRule(ASSAY_MAP_FILE_NAME, true, "The name of the file containing the assay to serum map"),
                AttributeRule.newIntegerRule(MDS_DIMENSION, false, "The dimension of the space for MDS"),
                AttributeRule.newBooleanRule(LOG_2_TRANSFORM, true, "Whether to log2 transform the data"),
                new ElementRule(TIP_TRAIT, CompoundParameter.class, "The parameter of tip locations from the tree", true),
                new ElementRule(LOCATIONS, MatrixParameter.class),
                new ElementRule(MDS_PRECISION, Parameter.class)
        };

        public Class getReturnType() {
            return AntigenicTraitLikelihood.class;
        }
    };

    public List<Citation> getCitations() {
        List<Citation> citations = new ArrayList<Citation>();
        citations.add(new Citation(
                new Author[]{
                        new Author("A", "Rambaut"),
                        new Author("T", "Bedford"),
                        new Author("P", "Lemey"),
                        new Author("C", "Russell"),
                        new Author("D", "Smith"),
                        new Author("MA", "Suchard"),
                },
                Citation.Status.IN_PREPARATION
        ));
        return citations;
    }
}

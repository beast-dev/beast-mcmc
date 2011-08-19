package dr.evomodel.continuous;

import dr.inference.model.*;
import dr.math.MathUtils;
import dr.util.*;
import dr.xml.*;

import java.io.*;
import java.util.*;
import java.util.logging.Logger;

/**
 * @author Andrew Rambaut
 * @author Marc Suchard
 * @version $Id$
 */
public class DiscreteAntigenicTraitLikelihood extends AntigenicTraitLikelihood implements Citable {

    public final static String DISCRETE_ANTIGENIC_TRAIT_LIKELIHOOD = "discreteAntigenicTraitLikelihood";

    /**
     * Constructor
     * @param mdsDimension dimension of the mds space
     * @param mdsPrecision parameter which gives the precision of the bmds
     * @param locationsParameter a parameter of locations of viruses/sera
     * @param dataTable the assay table (virus in rows, serum assays in columns)
     * @param virusAntiserumMap a map of viruses to corresponding sera
     * @param assayAntiserumMap a map of repeated assays for a given sera
     * @param log2Transform transform the data into log 2 space
     */
    public DiscreteAntigenicTraitLikelihood(
            int mdsDimension,
            Parameter mdsPrecision,
            Parameter clusterIndexParameter,
            MatrixParameter locationsParameter,
            DataTable<String[]> dataTable,
            Map<String, String> virusAntiserumMap,
            Map<String, String> assayAntiserumMap,
            final boolean log2Transform) {

        super(mdsDimension,
                mdsPrecision,
                null,
                locationsParameter,
                dataTable,
                virusAntiserumMap,
                assayAntiserumMap,
                log2Transform);

        // Start off with a 1-to-1 correspondence between location and cluster
//        maxClusterCount = getLocationCount();

        maxClusterCount = 40;

        this.clusterIndexParameter = clusterIndexParameter;

        clusterIndexParameter.setDimension(getLocationCount());
        clusterSizes = new int[maxClusterCount];

        //Force the boundaries of rateCategoryParameter to match the category count
        Parameter.DefaultBounds bound = new Parameter.DefaultBounds(maxClusterCount - 1, 0, maxClusterCount);
        clusterIndexParameter.addBounds(bound);

        for (int i = 0; i < maxClusterCount; i++) {
            int r = MathUtils.nextInt(maxClusterCount);
            clusterIndexParameter.setParameterValue(i, r);
        }

        updateClusterSizes();

        addVariable(clusterIndexParameter);
        addStatistic(new ClusterMask());
        addStatistic(new ClusterIndices());
        addStatistic(new ClusterCount());
        addStatistic(new ClusterSizes());
        addStatistic(new ClusteredLocations());
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Variable.ChangeType type) {
        if (variable == clusterIndexParameter) {
            if (index != -1) {
                // one dimension has changed
                int j = (int)clusterIndexParameter.getParameterValue(index);
                locationUpdated[j] = true;

            } else {
                // all have changed
                for (int i = 0; i < maxClusterCount; i++) {
                    int j = (int)clusterIndexParameter.getParameterValue(i);
                    locationUpdated[j] = true;
                }
            }
            distancesKnown = false;
            thresholdsKnown = false;
            truncationKnown = false;

            clusterMaskKnown = false;
        }

        super.handleVariableChangedEvent(variable, index, type);
    }

    @Override
    public void makeDirty() {
        super.makeDirty();
        clusterMaskKnown = false;
    }

    Double[] values;

    @Override
    protected void storeState() {
        super.storeState();

        values = clusterIndexParameter.getValues();
    }

    @Override
    protected void restoreState() {
        super.restoreState();

        clusterMaskKnown = false;
    }

    @Override
    protected int getLocationIndex(final int index) {
        return (int)clusterIndexParameter.getParameterValue(index);
    }

    private void updateClusterSizes() {
        for (int i = 0; i < maxClusterCount; i++) {
            clusterSizes[i] = 0;
        }
        for (int i = 0; i < maxClusterCount; i++) {
            int j = (int)clusterIndexParameter.getParameterValue(i);
            clusterSizes[j] ++;
        }
        clusterCount = 0;
        for (int i = 0; i < maxClusterCount; i++) {
            if (clusterSizes[i] > 0) {
                clusterCount++;
            }
        }
        clusterMaskKnown = true;
    }

    private int maxClusterCount;
    private final Parameter clusterIndexParameter;
    private final int[] clusterSizes;
    private int clusterCount;

    private boolean clusterMaskKnown;

    public class ClusterMask extends Statistic.Abstract {

        public ClusterMask() {
            super("clusterMask");
        }

        public int getDimension() {
            return maxClusterCount;
        }

        public double getStatisticValue(int i) {
            if (!clusterMaskKnown) {
                updateClusterSizes();
            }
            return clusterSizes[i] > 0 ? 1.0 : 0.0;
        }

    }

    public class ClusterIndices extends Statistic.Abstract {

        public ClusterIndices() {
            super("clusterIndices");
        }

        public int getDimension() {
            return clusterIndexParameter.getDimension();
        }

        @Override
        public String getDimensionName(final int i) {
            Parameter loc = getLocationsParameter().getParameter(i);
            return loc.getParameterName();
        }

        public double getStatisticValue(int i) {
            return clusterIndexParameter.getParameterValue(i);
        }

    }

    public class ClusterCount extends Statistic.Abstract {

        public ClusterCount() {
            super("clusterCount");
        }

        public int getDimension() {
            return 1;
        }

        public double getStatisticValue(int i) {
            if (!clusterMaskKnown) {
                updateClusterSizes();
            }
            return clusterCount;
        }

    }

    public class ClusterSizes extends Statistic.Abstract {

        public ClusterSizes() {
            super("clusterSizes");
        }

        public int getDimension() {
            return maxClusterCount;
        }

        public double getStatisticValue(int i) {
            if (!clusterMaskKnown) {
                updateClusterSizes();
            }
            return clusterSizes[i];
        }

    }

    public class ClusteredLocations extends Statistic.Abstract {

        public ClusteredLocations() {
            super("clusterLocations");
        }

        @Override
        public String getDimensionName(final int i) {
            int location = i / getMDSDimension();
            int dim = i % getMDSDimension();

            Parameter loc = getLocationsParameter().getParameter(location);
            if (getMDSDimension() == 2) {
                return loc.getParameterName() + "_" + (dim == 0 ? "X" : "Y");
            } else {
                return loc.getParameterName() + "_" + (dim + 1);
            }
        }

        public int getDimension() {
            return maxClusterCount * getMDSDimension();
        }

        public double getStatisticValue(final int i) {
            int location = i / getMDSDimension();
            int dim = i % getMDSDimension();

            int j = (int)clusterIndexParameter.getParameterValue(location);
            Parameter loc = getLocationsParameter().getParameter(j);

            return loc.getParameterValue(dim);
        }

    }

    // **************************************************************
    // XMLObjectParser
    // **************************************************************

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {
        public final static String FILE_NAME = "fileName";
        public final static String VIRUS_MAP_FILE_NAME = "virusMapFile";
        public final static String ASSAY_MAP_FILE_NAME = "assayMapFile";

        public static final String CLUSTER_INDICES = "clusterIndices";
        public final static String TIP_TRAIT = "tipTrait";
        public final static String LOCATIONS = "locations";
        public static final String MDS_DIMENSION = "mdsDimension";
        public static final String MDS_PRECISION = "mdsPrecision";

        public static final String LOG_2_TRANSFORM = "log2Transform";
        public static final String TITRATION_THRESHOLD = "titrationThreshold";

        public String getParserName() {
            return DISCRETE_ANTIGENIC_TRAIT_LIKELIHOOD;
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

            Parameter clusterIndicesParameter = (Parameter) xo.getElementFirstChild(CLUSTER_INDICES);

            // This parameter needs to be linked to the one in the IntegratedMultivariateTreeLikelihood (I suggest that the parameter is created
            // here and then a reference passed to IMTL - which optionally takes the parameter of tip trait values, in which case it listens and
            // updates accordingly.
            CompoundParameter tipTraitParameter = null;
            if (xo.hasChildNamed(TIP_TRAIT)) {
                tipTraitParameter = (CompoundParameter) xo.getElementFirstChild(TIP_TRAIT);
            }

            MatrixParameter locationsParameter = (MatrixParameter) xo.getElementFirstChild(LOCATIONS);

            Parameter mdsPrecision = (Parameter) xo.getElementFirstChild(MDS_PRECISION);

            AntigenicTraitLikelihood AGTL = new DiscreteAntigenicTraitLikelihood(mdsDimension, mdsPrecision, clusterIndicesParameter, locationsParameter, assayTable, virusAntiserumMap, assayAntiserumMap, log2Transform);

            Logger.getLogger("dr.evomodel").info("Using Discrete Evolutionary Cartography model. Please cite:\n" + Utils.getCitationString(AGTL));

            return AGTL;
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
                    "for viruses and sera/antisera in some multidimensional 'antigenic' space. This is a discrete classifier form of the model" +
                    "which assigns viruses to discrete antigenic classes.";
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
                new ElementRule(CLUSTER_INDICES, Parameter.class, "The parameter of cluster indices for each virus/serum"),
                new ElementRule(TIP_TRAIT, CompoundParameter.class, "The parameter of tip locations from the tree", true),
                new ElementRule(LOCATIONS, MatrixParameter.class),
                new ElementRule(MDS_PRECISION, Parameter.class)
        };

        public Class getReturnType() {
            return DiscreteAntigenicTraitLikelihood.class;
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

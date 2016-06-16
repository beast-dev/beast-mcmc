/*
 * DiscreteAntigenicTraitLikelihood.java
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
@Deprecated
public class DiscreteAntigenicTraitLikelihood extends AntigenicTraitLikelihood implements Citable {

    public final static String DISCRETE_ANTIGENIC_TRAIT_LIKELIHOOD = "discreteAntigenicTraitLikelihood";

    private static final int CLUSTER_COUNT = -1;

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
            CompoundParameter tipTraitParameter,
            DataTable<String[]> dataTable,
            Map<String, String> virusAntiserumMap,
            Map<String, String> assayAntiserumMap,
            List<String> virusLocationStatisticList,
            final boolean log2Transform) {

        super(DISCRETE_ANTIGENIC_TRAIT_LIKELIHOOD);

        String[] virusNames = dataTable.getRowLabels();
        String[] assayNames = dataTable.getColumnLabels();

        // the total number of viruses is the number of rows in the table
        int virusCount = dataTable.getRowCount();
        int assayCount = dataTable.getColumnCount();

        int[] assayToSerumIndices = new int[assayNames.length];

        double[][] observationValueTable = new double[virusCount][assayCount];
        ObservationType[][] observationTypeTable = new ObservationType[virusCount][assayCount];

        initalizeTable(dataTable, observationValueTable, observationTypeTable, log2Transform);

        // This removes viruses that are not in the tree
        List<String> tipLabels = null;
        if (tipTraitParameter != null) {
            tipLabels = new ArrayList<String>();
            int tipCount = tipTraitParameter.getParameterCount();

            for (int i = 0; i < tipCount; i++) {
                String label = tipTraitParameter.getParameter(i).getParameterName();
                if (label.endsWith(".antigenic")) {
                    label = label.substring(0, label.indexOf(".antigenic"));
                }
                tipLabels.add(label);
            }
        }

        // locations are either viruses or sera (or both)
        List<String> locationLabelsList = new ArrayList<String>();
        int[] virusToLocationIndices = new int[virusCount];
        int count = 0;
        for (String virusName : virusNames) {
            String name = null;

            if (virusAntiserumMap != null) {
                name = virusAntiserumMap.get(virusName);
            }

            if (name == null) {
                name = virusName;
            }

            virusToLocationIndices[count] = locationLabelsList.size();
            locationLabelsList.add(name);

            count++;
        }

        List<String> serumNamesList = new ArrayList<String>();
        count = 0;
        for (String assayName : assayNames) {
            String name = null;

            if (assayAntiserumMap != null) {
                name = assayAntiserumMap.get(assayName);
            }

            if (name == null) {
                name = assayName;
            }

            int index = serumNamesList.indexOf(name);
            if (index == -1) {
                index = serumNamesList.size();
                serumNamesList.add(name);
            }
            assayToSerumIndices[count] = index;
            count++;
        }
        String[] serumNames = new String[serumNamesList.size()];
        serumNamesList.toArray(serumNames);

        int serumCount = serumNames.length;
        int[] serumToLocationIndices = new int[serumCount];
        count = 0;
        for (String serumName : serumNames) {
            int index = locationLabelsList.indexOf(serumName);
            if (index == -1) {
                index = locationLabelsList.size();
                locationLabelsList.add(serumName);
            }
            serumToLocationIndices[count] = index;
            count++;
        }

        String[] locationLabels = new String[locationLabelsList.size()];
        locationLabelsList.toArray(locationLabels);
        int locationCount = locationLabels.length;

        List<Double> observationList = new ArrayList<Double>();
        List<ObservationType> observationTypeList = new ArrayList<ObservationType>();

        int[] virusObservationCounts = new int[virusCount];
        int[] serumObservationCounts = new int[serumCount];

        List<Pair> locationPairs = new ArrayList<Pair>();

//        System.out.print("virus");
//        for (int j = 0; j < assayCount; j++) {
//            int k = assayToSerumIndices[j];
//            System.out.print("\t" + assayNames[k] + "[" + locationLabels[serumToLocationIndices[k]] + "]");
//        }
//        System.out.println();

        // Build a sparse matrix of non-missing assay values
        for (int i = 0; i < virusCount; i++) {

            if (virusToLocationIndices[i] != -1) {
                // viruses with location indices of minus one have been excluded

//                System.out.print(virusNames[i] + "[" + locationLabels[virusToLocationIndices[i]] + "]");

                for (int j = 0; j < assayCount; j++) {
                    int k = assayToSerumIndices[j];

                    Double value = observationValueTable[i][j];
                    ObservationType type = observationTypeTable[i][j];

//                    System.out.print("\t" + value);

                    if (type != ObservationType.MISSING) {
                        observationList.add(value);
                        observationTypeList.add(type);

                        locationPairs.add(new Pair(virusToLocationIndices[i], serumToLocationIndices[k]));

                        virusObservationCounts[i]++;
                        serumObservationCounts[k]++;
                    }
                }
//                System.out.println();
            }
        }

        // check that all the viruses and sera have observations
        for (int i = 0; i < virusCount; i++) {
            if (virusToLocationIndices[i] != -1 && virusObservationCounts[i] == 0) {
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

        int[] rowLocationIndices = new int[locationPairs.size()];
        for (int i = 0; i < rowLocationIndices.length; i++) {
            rowLocationIndices[i] = locationPairs.get(i).location1;
        }

        int[] columnLocationIndices = new int[locationPairs.size()];
        for (int i = 0; i < columnLocationIndices.length; i++) {
            columnLocationIndices[i] = locationPairs.get(i).location2;
        }

        ObservationType[] observationTypes = new ObservationType[observationTypeList.size()];
        observationTypeList.toArray(observationTypes);

        int thresholdCount = 0;
        for (int i = 0; i < observations.length; i++) {
            thresholdCount += (observationTypes[i] != ObservationType.POINT ? 1 : 0);
        }

        if (tipTraitParameter != null) {
            //  the location -> tip map
            tipIndices = new int[locationCount];

            for (int i = 0; i < locationCount; i++) {
                tipIndices[i] = tipLabels.indexOf(locationLabels[i]);
            }

//            for (int i = 0; i < locationCount; i++) {
//                if (tipIndices[i] == -1) {
//                    System.err.println("Location, " + locationLabels[i] + ", not found in tree");
//                }
//            }

            for (String tipLabel : tipLabels) {
                if (!locationLabelsList.contains(tipLabel)) {
                    System.err.println("Tip, " + tipLabel + ", not found in location list");
                }
            }
        } else {
            tipIndices = null;
        }

        this.tipTraitParameter = tipTraitParameter;

        StringBuilder sb = new StringBuilder();
        sb.append("\tDiscreteAntigenicTraitLikelihood:\n");
        sb.append("\t\t" + virusNames.length + " viruses\n");
        sb.append("\t\t" + assayNames.length + " assays\n");
        sb.append("\t\t" + serumNames.length + " antisera\n");
        sb.append("\t\t" + locationLabels.length + " locations\n");
        sb.append("\t\t" + locationPairs.size() + " distances\n");
        sb.append("\t\t" + observations.length + " observations\n");
        sb.append("\t\t" + thresholdCount + " threshold observations\n");
        Logger.getLogger("dr.evomodel").info(sb.toString());

        initialize(
                mdsDimension,
                false,
                mdsPrecision,
                locationsParameter,
                locationLabels,
                observations,
                observationTypes,
                rowLocationIndices,
                columnLocationIndices);

        // some random initial locations
        for (int i = 0; i < locationsParameter.getParameterCount(); i++) {
            for (int j = 0; j < mdsDimension; j++) {
             //   double r = MathUtils.nextGaussian();
                double r = 0.0;
                if (j == 0) {
                    r = (double) i * 0.05;
                }
                else {
                    r = MathUtils.nextGaussian();
                }
                locationsParameter.getParameter(i).setParameterValueQuietly(j, r);
            }
        }

        // Start off with a 1-to-1 correspondence between location and cluster
        if (CLUSTER_COUNT > 0) {
            maxClusterCount = CLUSTER_COUNT;
        } else {
            maxClusterCount = getLocationCount();
        }

        this.clusterIndexParameter = clusterIndexParameter;

        clusterIndexParameter.setDimension(getLocationCount());
        clusterSizes = new int[maxClusterCount];

        //Force the boundaries of rateCategoryParameter to match the category count
        Parameter.DefaultBounds bound = new Parameter.DefaultBounds(maxClusterCount - 1, 0, getLocationCount());
        clusterIndexParameter.addBounds(bound);

        for (int i = 0; i < getLocationCount(); i++) {
            int r = i;
            clusterIndexParameter.setParameterValue(i, r);

            for (int dim = 0; dim < mdsDimension; dim++) {
                tipTraitParameter.setParameterValue((tipIndices[r] * mdsDimension) + dim,
                        locationsParameter.getParameterValue((r * mdsDimension) + dim));
            }
        }

        updateClusterSizes();

        addVariable(clusterIndexParameter);
        addStatistic(new ClusterMask());
        addStatistic(new ClusterIndices());
        addStatistic(new ClusterCount());
        addStatistic(new ClusterSizes());
        addStatistic(new ClusteredLocations());

        int i = 0;
        for (String virusName : virusNames) {
            if (virusLocationStatisticList.contains(virusName)) {
                addStatistic(new VirusLocation(virusName + "." + "location", i));
            }
            i++;
        }
    }

    @Override
    protected void setupLocationsParameter(MatrixParameter locationsParameter) {
        locationsParameter.setColumnDimension(getMDSDimension());
        int n = CLUSTER_COUNT;
        if (n < 1) {
            n = getLocationCount();
        }
        locationsParameter.setRowDimension(n);
        for (int i = 0; i < n; i++) {
            locationsParameter.getParameter(i).setId("cluster_" + (i+1));
        }
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Variable.ChangeType type) {
        if (variable == clusterIndexParameter) {
            for (int i = 0; i < distanceUpdated.length; i++) {
                distanceUpdated[i] = true;
            }
            residualsKnown = false;
            thresholdsKnown = false;
            clusterMaskKnown = false;
        }

        if (tipTraitParameter != null) {
            Parameter locations = getLocationsParameter();
            int mdsDimension = getMDSDimension();

            if  (variable == locations) {
                int location = index / mdsDimension;
                int dim = index % mdsDimension;

                for (int i = 0; i < clusterIndexParameter.getDimension(); i++) {
                    // the location of a cluster has moved so set the new location of every virus in that cluster
                    if (((int)clusterIndexParameter.getParameterValue(i)) == location) {
                        if (tipIndices[i] != -1) {
                            tipTraitParameter.setParameterValue((tipIndices[i] * mdsDimension) + dim, locations.getParameterValue(index));
                        }
                    }
                }
            } else if (variable == clusterIndexParameter) {
                // a virus has moved cluster so set its location to that of a new cluster
                if (tipIndices[index] != -1) {
                    int location = (int)clusterIndexParameter.getParameterValue(index);
                    for (int dim = 0; dim < mdsDimension; dim++) {
                        tipTraitParameter.setParameterValue((tipIndices[index] * mdsDimension) + dim, locations.getParameterValue((location * mdsDimension) + dim));
                    }
                }
            }
        }

        super.handleVariableChangedEvent(variable, index, type);
    }

    public CompoundParameter getTipTraitParameter() {
        return tipTraitParameter;
    }

    public int[] getTipIndices() {
        return tipIndices;
    }

    private CompoundParameter tipTraitParameter;
    private int[] tipIndices;

    @Override
    public void makeDirty() {
        super.makeDirty();
        clusterMaskKnown = false;
    }

    @Override
    protected void storeState() {
        super.storeState();
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
        for (int i = 0; i < getLocationCount(); i++) {
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
            return getLocationLabels()[i];
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
            super("clusteredLocations");
        }

        @Override
        public String getDimensionName(final int i) {
            int location = i / getMDSDimension();
            int dim = i % getMDSDimension();

            String label = getLocationLabels()[location];
            if (getMDSDimension() == 2) {
                return label + "_" + (dim == 0 ? "X" : "Y");
            } else {
                return label + "_" + (dim + 1);
            }
        }

        public int getDimension() {
            return getLocationCount() * getMDSDimension();
        }

        public double getStatisticValue(final int i) {
            int location = i / getMDSDimension();
            int dim = i % getMDSDimension();

            int j = (int)clusterIndexParameter.getParameterValue(location);
            Parameter loc = getLocationsParameter().getParameter(j);

            return loc.getParameterValue(dim);
        }

    }

    public class VirusLocation extends Statistic.Abstract {

        public VirusLocation(String statisticName, int virusIndex) {
            super(statisticName);
            this.virusIndex = virusIndex;
        }

        @Override
        public String getDimensionName(final int dim) {
            if (getMDSDimension() == 2) {
                return getStatisticName() + "_" + (dim == 0 ? "X" : "Y");
            } else {
                return getStatisticName() + "_" + (dim + 1);
            }
        }

        public int getDimension() {
            return getMDSDimension();
        }

        public double getStatisticValue(final int dim) {
            int cluster = (int)clusterIndexParameter.getParameterValue(virusIndex);
            Parameter loc = getLocationsParameter().getParameter(cluster);
            return loc.getParameterValue(dim);
        }

        private final int virusIndex;
    }

    private class Pair {
        Pair(final int location1, final int location2) {
            if (location1 < location2) {
                this.location1 = location1;
                this.location2 = location2;
            } else {
                this.location1 = location2;
                this.location2 = location1;
            }
        }

        int location1;
        int location2;

        @Override
        public boolean equals(final Object o) {
            return ((Pair)o).location1 == location1 && ((Pair)o).location2 == location2;
        }

        @Override
        public String toString() {
            return "" + location1 + ", " + location2;
        }
    };
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

        public static final String VIRUS_LOCATIONS = "virusLocations";

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
                throw new XMLParseException("Unable to read assay data from file: " + e.getMessage());
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

            List<String> virusLocationStatisticList = null;
            String[] virusLocations = xo.getStringArrayAttribute(VIRUS_LOCATIONS);
            if (virusLocations != null) {
                virusLocationStatisticList = Arrays.asList(virusLocations);
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

            AntigenicTraitLikelihood AGTL = new DiscreteAntigenicTraitLikelihood(mdsDimension, mdsPrecision, clusterIndicesParameter, locationsParameter, tipTraitParameter, assayTable, virusAntiserumMap, assayAntiserumMap, virusLocationStatisticList, log2Transform);

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
                AttributeRule.newStringArrayRule(VIRUS_LOCATIONS, true, "A list of virus names to create location statistics for"),
                new ElementRule(CLUSTER_INDICES, Parameter.class, "The parameter of cluster indices for each virus/serum"),
                new ElementRule(TIP_TRAIT, CompoundParameter.class, "The parameter of tip locations from the tree", true),
                new ElementRule(LOCATIONS, MatrixParameter.class),
                new ElementRule(MDS_PRECISION, Parameter.class)
        };

        public Class getReturnType() {
            return DiscreteAntigenicTraitLikelihood.class;
        }
    };

    @Override
    public Citation.Category getCategory() {
        return Citation.Category.TRAIT_MODELS;
    }

    @Override
    public String getDescription() {
        return "Bayesian Antigenic Cartography framework";
    }

    public List<Citation> getCitations() {
        return Arrays.asList(CommonCitations.BEDFORD_2015_INTEGRATING);
    }
}

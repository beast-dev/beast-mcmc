/*
 * ContinuousAntigenicTraitLikelihood.java
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
 * Antigenic evolution assuming each virus's antigenic property diffuses independently.
 * @author Andrew Rambaut
 * @author Marc Suchard
 * @version $Id$
 */
@Deprecated // for the moment at least
public class ContinuousAntigenicTraitLikelihood extends AntigenicTraitLikelihood implements Citable {

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
    public ContinuousAntigenicTraitLikelihood(
            int mdsDimension,
            Parameter mdsPrecision,
            CompoundParameter tipTraitParameter,
            MatrixParameter locationsParameter,
            DataTable<String[]> dataTable,
            Map<String, String> virusAntiserumMap,
            Map<String, String> assayAntiserumMap,
            List<String> virusLocationStatisticList,
            final boolean log2Transform) {

        super(ANTIGENIC_TRAIT_LIKELIHOOD);

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

            if (tipLabels == null || tipLabels.contains(name)) {
                virusToLocationIndices[count] = locationLabelsList.size();
                locationLabelsList.add(name);
            } else {
                virusToLocationIndices[count] = -1;
            }
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

        // Build a sparse matrix of non-missing assay values
        for (int i = 0; i < virusCount; i++) {

            if (virusToLocationIndices[i] != -1) {
                // viruses with location indices of minus one have been excluded

                for (int j = 0; j < assayCount; j++) {
                    int k = assayToSerumIndices[j];

                    Double value = observationValueTable[i][j];
                    ObservationType type = observationTypeTable[i][j];

                    if (type != ObservationType.MISSING) {
                        observationList.add(value);
                        observationTypeList.add(type);

                        locationPairs.add(new Pair(virusToLocationIndices[i], serumToLocationIndices[k]));

                        virusObservationCounts[i]++;
                        serumObservationCounts[k]++;
                    }
                }
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

            for (int i = 0; i < locationCount; i++) {
                if (tipIndices[i] == -1) {
                    System.err.println("Location, " + locationLabels[i] + ", not found in tree");
                }
            }

            for (String tipLabel : tipLabels) {
                if (!locationLabelsList.contains(tipLabel)) {
                    System.err.println("Tip, " + tipLabel + ", not found in location list");
                }
            }
        } else {
            tipIndices = null;
        }

        StringBuilder sb = new StringBuilder();
        sb.append("\tAntigenicTraitLikelihood:\n");
        sb.append("\t\t" + virusNames.length + " viruses\n");
        sb.append("\t\t" + assayNames.length + " assays\n");
        sb.append("\t\t" + serumNames.length + " antisera\n");
        sb.append("\t\t" + locationLabels.length + " locations\n");
        sb.append("\t\t" + locationPairs.size() + " distances\n");
        sb.append("\t\t" + observations.length + " observations\n");
        sb.append("\t\t" + thresholdCount + " threshold observations\n");
        Logger.getLogger("dr.evomodel").info(sb.toString());

        this.tipTraitParameter = tipTraitParameter;

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

                if (tipTraitParameter != null) {
                    if (tipIndices[i] != -1) {
                        tipTraitParameter.setParameterValue((tipIndices[i] * mdsDimension) + j, r);
                    }
                }
            }
        }

        int i = 0;
        for (String virusName : virusNames) {
	        if (virusLocationStatisticList.contains(virusName)) {
                addStatistic(new VirusLocation(virusName + "." + "location", i));
	        }
            i++;
        }

    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Variable.ChangeType type) {
        if (tipTraitParameter != null && variable == getLocationsParameter()) {
            int mdsDimension = getMDSDimension();

            int location = index / mdsDimension;
            int dim = index % mdsDimension;

            if (tipIndices[location] != -1) {
                double value = getLocationsParameter().getParameterValue(index);
                tipTraitParameter.setParameterValue((tipIndices[location] * mdsDimension) + dim, value);
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
            Parameter loc = getLocationsParameter().getParameter(virusIndex);
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
    };

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

        public static final String VIRUS_LOCATIONS = "virusLocations";

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

            ContinuousAntigenicTraitLikelihood AGTL = new ContinuousAntigenicTraitLikelihood(mdsDimension, mdsPrecision, tipTraitParameter, locationsParameter, assayTable, virusAntiserumMap, assayAntiserumMap, virusLocationStatisticList, log2Transform);

            Logger.getLogger("dr.evomodel").info("Using EvolutionaryCartography model. Please cite:\n" + Utils.getCitationString(AGTL));

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
                AttributeRule.newStringArrayRule(VIRUS_LOCATIONS, true, "A list of virus names to create location statistics for"),
                new ElementRule(TIP_TRAIT, CompoundParameter.class, "The parameter of tip locations from the tree", true),
                new ElementRule(LOCATIONS, MatrixParameter.class),
                new ElementRule(MDS_PRECISION, Parameter.class)
        };

        public Class getReturnType() {
            return ContinuousAntigenicTraitLikelihood.class;
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

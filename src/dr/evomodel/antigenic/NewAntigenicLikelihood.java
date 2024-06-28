/*
 * NewAntigenicLikelihood.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.evomodel.antigenic;

import dr.inference.model.*;
import dr.inference.multidimensionalscaling.*;
import dr.math.MathUtils;
import dr.util.Citable;
import dr.util.Citation;
import dr.util.CommonCitations;
import dr.util.DataTable;

import java.util.*;
import java.util.logging.Logger;


/**
 * @author Andrew Rambaut
 * @author Trevor Bedford
 * @author Marc Suchard
 */
/*
    Both virus locations and serum locations are shifted by the parameter locationDrift.
    A location is increased by locationDrift x offset.
    Offset is set to 0 for the earliest virus and increasing with difference in date from the earliest virus.
    This makes the raw virusLocations and serumLocations parameters not directly interpretable.
*/
public class NewAntigenicLikelihood extends AbstractModelLikelihood implements Citable {

//    private static final boolean CHECK_INFINITE = false;
    private static final boolean USE_THRESHOLDS = true;
    private static final boolean USE_INTERVALS = true;

    public final static String ANTIGENIC_LIKELIHOOD = "newAntigenicLikelihood";

    // column indices in table
//    private static final int VIRUS_ISOLATE = 0;
    private static final int VIRUS_STRAIN = 1;
    private static final int VIRUS_DATE = 2;
    private static final int SERUM_ISOLATE = 3;
    private static final int SERUM_STRAIN = 4;
    private static final int SERUM_DATE = 5;
    private static final int TITRE = 6;

    public enum MeasurementType {
        INTERVAL,
        POINT,
        THRESHOLD,
        MISSING
    }

    public NewAntigenicLikelihood(
            int mdsDimension,
            Parameter mdsPrecisionParameter,
            Parameter locationDriftParameter,
            Parameter virusDriftParameter,
            Parameter serumDriftParameter,
            MatrixParameter virusSamplingParameter, // TODO Remove
            MatrixParameterInterface serumLocationsParameter,
            CompoundParameter tipTraitsParameter,
            Parameter virusOffsetsParameter,
            Parameter serumOffsetsParameter,
            Parameter serumPotenciesParameter,
            Parameter serumBreadthsParameter,
            Parameter virusAviditiesParameter,
            DataTable<String[]> dataTable,
            boolean mergeSerumIsolates,
            double intervalWidth,
            double driftInitialLocations,
            int tipStartOffset) {

        super(ANTIGENIC_LIKELIHOOD);

        this.intervalWidth = intervalWidth;
        boolean useIntervals = USE_INTERVALS && intervalWidth > 0.0;

        MatchInfo info = matchVirusesAndSerum(measurements,
                dataTable,
                virusNames,
                virusDates,
                serumNames,
                serumDates,
                mergeSerumIsolates,
                useIntervals);

        double[] maxColumnTitres = getMaxTitres(measurements, serumNames);

        this.tipTraitsParameter = tipTraitsParameter;
        if (tipTraitsParameter.getBounds() == null) {
            tipTraitsParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY,
                    Double.NEGATIVE_INFINITY, tipTraitsParameter.getDimension()));
        }
        addVariable(tipTraitsParameter);
        this.tipIndices = setupTipIndices(this.tipTraitsParameter, virusNames);
//        this.virusIndices = setupVirusIndices(tipIndices);

        this.mdsDimension = mdsDimension;
        this.tipDimension = tipTraitsParameter.getParameter(0).getDimension();
        this.tipStartOffset = tipStartOffset;

        this.mdsPrecisionParameter = mdsPrecisionParameter;
        addVariable(mdsPrecisionParameter);

        this.locationDriftParameter = locationDriftParameter;
        if (this.locationDriftParameter != null) {
            addVariable(locationDriftParameter);

            if (virusOffsetsParameter == null || serumOffsetsParameter == null) {
                throw new IllegalArgumentException("Must also provide parameters to hold the offsets");
            }
        }

        this.virusDriftParameter = virusDriftParameter;
        if (this.virusDriftParameter != null) {
            addVariable(virusDriftParameter);
        }

        this.serumDriftParameter = serumDriftParameter;
        if (this.serumDriftParameter != null) {
            addVariable(serumDriftParameter);
        }

        this.virusSamplingParameter = virusSamplingParameter;
        if (this.virusSamplingParameter != null) {
            setupLocationsParameter(virusSamplingParameter, virusNames);
        }

        this.serumLocationsParameter = serumLocationsParameter;
        if (this.serumLocationsParameter != null) {
            setupLocationsParameter(serumLocationsParameter, serumNames);
        }

        this.virusOffsetsParameter = virusOffsetsParameter;
        if (virusOffsetsParameter != null) {
            setupOffsetsParameter(virusOffsetsParameter, virusNames, virusDates, info.earliestDate, true);
        }

        this.serumOffsetsParameter = serumOffsetsParameter;
        if (serumOffsetsParameter != null) {
            setupOffsetsParameter(serumOffsetsParameter, serumNames, serumDates, info.earliestDate, false);
        }

        this.serumPotenciesParameter = setupSerumPotencies(serumPotenciesParameter, maxColumnTitres);
        this.serumBreadthsParameter = setupSerumBreadths(serumBreadthsParameter);
        this.virusAviditiesParameter = setupVirusAvidities(virusAviditiesParameter);

        StringBuilder sb = new StringBuilder();
        sb.append("\tNewAntigenicLikelihood:\n");
        sb.append("\t\t").append(virusNames.size()).append(" viruses\n");
        sb.append("\t\t").append(serumNames.size()).append(" sera\n");
        sb.append("\t\t").append(measurements.size()).append(" assay measurements\n");
        if (USE_THRESHOLDS) {
            sb.append("\t\t").append(info.thresholdCount).append(" thresholded measurements\n");
        }
        if (useIntervals) {
            sb.append("\n\t\tAssuming a log 2 measurement interval width of ").append(intervalWidth).append("\n");
        }
        Logger.getLogger("dr.evomodel").info(sb.toString());

        numViruses = virusNames.size();
        numSera = serumNames.size();
        layout = (numViruses >= numSera) ?
                new Layout.VirusXSerum(numViruses, numSera, mdsDimension) :
                new Layout.SerumXVirus(numViruses, numSera, mdsDimension);

        layout.sort(measurements, tipIndices);

//        virusLocationChanged = new boolean[this.virusSamplingParameter.getParameterCount()];
//        serumLocationChanged = new boolean[this.serumLocationsParameter.getParameterCount()];
//        virusEffectChanged = new boolean[virusNames.size()];
//        serumEffectChanged = new boolean[serumNames.size()];
//        logLikelihoods = new double[measurements.size()];
//        storedLogLikelihoods = new double[measurements.size()];

        setupInitialLocations(driftInitialLocations);

        mdsCore = instantiateCore(mdsDimension, layout.getMajorDim(), layout.getMinorDim(), false);
        internalDimension = mdsCore.getInternalDimension();

        observationCount = transferObservations();
        mdsCore.setNonMissingObservationCount(observationCount);

        transferLocations();
        transferPrecision();

        makeDirty();
    }

    public int getNumberOfViruses() { return numViruses; }

    public int getNumberOfSera() { return numSera; }

    public int getMdsDimension() { return mdsDimension; }

    interface Layout {

        int getMajorDim();

        int getMinorDim();

        int getObservationIndex(int virus, int serum);

        int getVirusLocationOffset();

        int getSerumLocationOffset();

        void sort(List<Measurement> measurements, int[] virusIndices);

        abstract class Base implements Layout {

            final int numViruses;
            final int numSera;
            final int mdsDim;

            Base(int numViruses, int numSera, int mdsDim) {
                this.numViruses = numViruses;
                this.numSera = numSera;
                this.mdsDim = mdsDim;
            }

            public void sort(List<Measurement> measurements, int[] virusIndices) {
                measurements.sort(getComparator(virusIndices));
            }

            abstract Comparator<Measurement> getComparator(int[] virusIndices);
        }

        class VirusXSerum extends Base {
            VirusXSerum(int numViruses, int numSera, int mdsDim) {
                super(numViruses, numSera, mdsDim);
            }

            @Override
            public int getMajorDim() {
                return numViruses;
            }

            @Override
            public int getMinorDim() {
                return numSera;
            }

            @Override
            public int getObservationIndex(int virus, int serum) {
                return virus * numSera + serum;
            }

            @Override
            public int getVirusLocationOffset() {
                return 0;
            }

            @Override
            public int getSerumLocationOffset() {
                return numViruses * mdsDim;
            }

            @Override
            Comparator<Measurement> getComparator(int[] virusIndices) {
                return (lhs, rhs) -> {
                    int virusDiff = virusIndices[lhs.virus] - virusIndices[rhs.virus];
                    if (virusDiff != 0) {
                        return virusDiff;
                    } else {
                        return lhs.serum - rhs.serum;
                    }
                };
            }
        }

        class SerumXVirus extends Base {
            SerumXVirus(int numViruses, int numSera, int mdsDim) {
                super(numViruses, numSera, mdsDim);
            }

            @Override
            public int getMajorDim() {
                return numSera;
            }

            @Override
            public int getMinorDim() {
                return numViruses;
            }

            @Override
            public int getObservationIndex(int virus, int serum) {
                return serum * numViruses + virus;
            }

            @Override
            public int getVirusLocationOffset() {
                return numSera * mdsDim;
            }

            @Override
            public int getSerumLocationOffset() {
                return 0;
            }

            @Override
            Comparator<Measurement> getComparator(int[] virusIndices) {
                return (lhs, rhs) -> {
                    int serumDiff = lhs.serum - rhs.serum;
                    if (serumDiff != 0) {
                        return serumDiff;
                    } else {
                        return virusIndices[lhs.virus] - virusIndices[rhs.virus];
                    }
                };
            }
        }
    }

    private final int numViruses;
    private final int numSera;

    private final int observationCount;
    private final Layout layout;

    private boolean locationsKnown;
    private boolean observationsKnown;
    private boolean internalGradientKnown;
    private boolean precisionKnown;

    private boolean savedLikelihoodKnown;
    private double savedLogLikelihood;
    private boolean savedLocationsKnown;
    private boolean savedObservationsKnown;
    private boolean savedPrecisionKnown;

    public void updateParametersOnDevice() {
        if (!precisionKnown) {
            transferPrecision();
            precisionKnown = true;
        }

        if (!locationsKnown) {
            transferLocations();
            locationsKnown = true;
        }

        if (!observationsKnown) {
            transferObservations();
            observationsKnown = true;
        }
    }

    private double computeLogLikelihoodMds() {

        updateParametersOnDevice();
        mdsCore.makeDirty();
        return mdsCore.calculateLogLikelihood();
    }

    private double getAdjustedValue(Measurement measurement) {
        return  calculateBaseline(measurement.virus, measurement.serum) - measurement.log2Titre;
    }

    private void transferPrecision() {
        double[] precision = new double[]{ mdsPrecisionParameter.getParameterValue(0) };
        mdsCore.setParameters(precision);
    }

    private int transferObservations() {

        double[] observations = new double[numViruses * numSera];
        Arrays.fill(observations, Double.NaN);

        int count = 0;
        for (Measurement measurement : measurements) {
            final int idx1 = tipIndices[measurement.virus];
            final int idx2 = measurement.serum;

            if (measurement.type != MeasurementType.POINT) {
                throw new RuntimeException("Not yet implemented");
            }

            observations[layout.getObservationIndex(idx1, idx2)] = getAdjustedValue(measurement);
            ++count;
        }

        mdsCore.setPairwiseData(observations);

        return count;
    }

    private void transferLocations() {

//        System.err.println(new WrappedVector.Raw(virusSamplingParameter.getParameterValues()));
//        System.err.println(new WrappedVector.Raw(tipTraitsParameter.getParameterValues()));

        double[] locations = new double[(numViruses + numSera) * mdsDimension];

        int offset = layout.getVirusLocationOffset();
        for (int i = 0; i < numViruses; ++i) {
            Parameter parameter = tipTraitsParameter.getParameter(i);
            for (int j = 0; j < mdsDimension; ++j) {
                locations[offset + j] = parameter.getParameterValue(tipStartOffset + j);
            }

            if (locationDriftParameter != null) {
                locations[offset] += locationDriftParameter.getParameterValue(0) *
                        virusOffsetsParameter.getParameterValue(i);
            }

            offset += mdsDimension;
        }

//        double[] debug = new double[numViruses * mdsDimension];
//        System.arraycopy(locations, layout.getVirusLocationOffset(), debug, 0, numViruses * mdsDimension);
//        System.err.println(new WrappedVector.Raw(debug));

        offset = layout.getSerumLocationOffset();
        for (int i = 0; i < numSera; ++i) {
            Parameter parameter = serumLocationsParameter.getParameter(i);
            for (int j = 0; j < mdsDimension; ++j) {
                locations[offset + j] = parameter.getParameterValue(j);
            }

            if (locationDriftParameter != null) {
                locations[offset] += locationDriftParameter.getParameterValue(0) *
                        serumOffsetsParameter.getParameterValue(i);
            }
            offset += mdsDimension;
        }

        mdsCore.updateLocation(-1, locations);
    }

    private final MultiDimensionalScalingCore mdsCore;
    private final int internalDimension;
//    private final long flags = 0L;

    public MultiDimensionalScalingCore getCore() {
        return mdsCore;
    }

    private MultiDimensionalScalingCore instantiateCore(int dim,
                                                        int rowLocationCount, int columnLocationCount,
                                                        boolean isLeftTruncated) {

        long flags = 0;
        String r = System.getProperty(MultiDimensionalScalingLikelihood.REQUIRED_FLAGS_PROPERTY);
        if (r != null) {
            flags = Long.parseLong(r.trim());
        }

        MultiDimensionalScalingCore core;
        if (flags >= MultiDimensionalScalingCore.USE_NATIVE_MDS) {
            System.err.println("Attempting to use a native MDS core with flag: " + flags + "; may the force be with you ....");
            core = new MassivelyParallelMDSImpl();
        } else {
            System.err.println("Compute mode found: " + flags);
            core = new MultiDimensionalScalingCoreImpl();
        }

        if (isLeftTruncated) {
            flags |= MultiDimensionalScalingCore.LEFT_TRUNCATION;
        }

        System.err.println("Initializing with flags: " + flags);

        core.initialize(dim, new MultiDimensionalScalingLayout(rowLocationCount, columnLocationCount), flags);

        return core;
    }

    private static double[] getMaxTitres(List<Measurement> measurements,
                                         List<String> serumNames) {
        double[] maxColumnTitres = new double[serumNames.size()];
        for (Measurement measurement : measurements) {
            double titre = measurement.log2Titre;
            if (Double.isNaN(titre)) {
                titre = measurement.log2Titre;
            }
            if (titre > maxColumnTitres[measurement.serum]) {
                maxColumnTitres[measurement.serum] = titre;
            }
        }
        return maxColumnTitres;
    }

    private static class MatchInfo {
        double earliestDate;
        int thresholdCount;

        MatchInfo(double earliestDate, int thresholdCount) {
            this.earliestDate = earliestDate;
            this.thresholdCount = thresholdCount;
        }
    }

    private static MatchInfo matchVirusesAndSerum(List<Measurement> measurements,
                                                  DataTable<String[]> dataTable,
                                                  List<String> virusNames,
                                                  List<Double> virusDates,
                                                  List<String> serumNames,
                                                  List<Double> serumDates,
                                                  boolean mergeSerumIsolates,
                                                  boolean useIntervals) {

        int thresholdCount = 0;

        double earliestDate = Double.POSITIVE_INFINITY;
        for (int i = 0; i < dataTable.getRowCount(); i++) {

            String[] values = dataTable.getRow(i);

            String virusName = values[VIRUS_STRAIN];
            double virusDate = Double.parseDouble(values[VIRUS_DATE]);
            int virus = virusNames.indexOf(virusName);
            if (virus == -1) {
                virusNames.add(virusName);
                virusDates.add(virusDate);
                virus = virusNames.size() - 1;
            }

            String serumName;
            if (mergeSerumIsolates) {
                serumName = values[SERUM_STRAIN];
            } else {
                serumName = values[SERUM_ISOLATE];
            }
            double serumDate = Double.parseDouble(values[SERUM_DATE]);
            int serum = serumNames.indexOf(serumName);
            if (serum == -1) {
                serumNames.add(serumName);
                serumDates.add(serumDate);
                serum = serumNames.size() - 1;
            }

            boolean isThreshold = false;
            boolean isLowerThreshold = false;
            double rawTitre = Double.NaN;
            if (values[TITRE].length() > 0) {
                try {
                    rawTitre = Double.parseDouble(values[TITRE]);
                } catch (NumberFormatException nfe) {
                    // check if threshold below
                    if (values[TITRE].contains("<")) {
                        rawTitre = Double.parseDouble(values[TITRE].replace("<",""));
                        isThreshold = true;
                        isLowerThreshold = true;
                        thresholdCount++;
                    }
                    // check if threshold above
                    if (values[TITRE].contains(">")) {
                        rawTitre = Double.parseDouble(values[TITRE].replace(">",""));
                        isThreshold = true;
                        isLowerThreshold = false;
                        thresholdCount++;
                        //throw new IllegalArgumentException("Error in measurement: unsupported greater than threshold at row " + (i+1));
                    }
                }
            }

            if (serumDate < earliestDate) {
                earliestDate = serumDate;
            }

            if (virusDate < earliestDate) {
                earliestDate = virusDate;
            }

            MeasurementType type = (isThreshold ? MeasurementType.THRESHOLD : (useIntervals ? MeasurementType.INTERVAL : MeasurementType.POINT));
            Measurement measurement = new Measurement(virus, serum, virusDate, serumDate, type, rawTitre, isLowerThreshold);

            if (USE_THRESHOLDS || !isThreshold) {
                measurements.add(measurement);
            }
        }

        return new MatchInfo(earliestDate, thresholdCount);
    }


    private Parameter setupVirusAvidities(Parameter virusAviditiesParameter) {
        // If no row parameter is given, then we will only use the serum effects
        if (virusAviditiesParameter != null) {
            virusAviditiesParameter.addBounds(new Parameter.DefaultBounds(Double.MAX_VALUE, Double.MIN_VALUE, 1));
            virusAviditiesParameter.setDimension(virusNames.size());
            addVariable(virusAviditiesParameter);
            String[] labelArray = new String[virusNames.size()];
            virusNames.toArray(labelArray);
            virusAviditiesParameter.setDimensionNames(labelArray);
            for (int i = 0; i < virusNames.size(); i++) {
                virusAviditiesParameter.setParameterValueQuietly(i, 0.0);
            }
        }
        return virusAviditiesParameter;
    }

    private Parameter setupSerumPotencies(Parameter serumPotenciesParameter, double[] maxColumnTitres) {
        // If no serum potencies parameter is given, make one to hold maximum values for scaling titres...
        if (serumPotenciesParameter == null) {
            serumPotenciesParameter = new Parameter.Default("serumPotencies");
        } else {
            serumPotenciesParameter.addBounds(new Parameter.DefaultBounds(Double.MAX_VALUE, 0.0, 1));
            addVariable(serumPotenciesParameter);
        }

        serumPotenciesParameter.setDimension(serumNames.size());
        String[] labelArray = new String[serumNames.size()];
        serumNames.toArray(labelArray);
        serumPotenciesParameter.setDimensionNames(labelArray);
        for (int i = 0; i < maxColumnTitres.length; i++) {
            serumPotenciesParameter.setParameterValueQuietly(i, maxColumnTitres[i]);
        }

        return serumPotenciesParameter;
    }

    private Parameter setupSerumBreadths(Parameter serumBreadthsParameter) {
        // If no serum breadths parameter is given, then we will only use the serum potencies
        if (serumBreadthsParameter != null) {
            serumBreadthsParameter.addBounds(new Parameter.DefaultBounds(Double.MAX_VALUE, 0.0, 1));
            serumBreadthsParameter.setDimension(serumNames.size());
            addVariable(serumBreadthsParameter);
            String[] labelArray = new String[serumNames.size()];
            serumNames.toArray(labelArray);
            serumBreadthsParameter.setDimensionNames(labelArray);
            for (int i = 0; i < serumNames.size(); i++) {
                serumBreadthsParameter.setParameterValueQuietly(i, 1.0);
            }
        }
        return serumBreadthsParameter;
    }

    protected void setupLocationsParameter(MatrixParameterInterface locationsParameter, List<String> strains) {
        if (locationsParameter instanceof MatrixParameter) {
            ((MatrixParameter) locationsParameter).setColumnDimension(mdsDimension);
            ((MatrixParameter) locationsParameter).setRowDimension(strains.size());
        } else if (locationsParameter instanceof FastMatrixParameter) {
            FastMatrixParameter fmp = (FastMatrixParameter) locationsParameter;
            if (fmp.getRowDimension() != mdsDimension) {
                throw new IllegalArgumentException("Column dim must be " + mdsDimension);
            }
            if (fmp.getColumnDimension() != strains.size()) {
                throw new IllegalArgumentException("Row dim must be " + strains.size());
            }
        }

        for (int i = 0; i < strains.size(); i++) {
            locationsParameter.getParameter(i).setId(strains.get(i));
        }
        addVariable(locationsParameter);
    }

    private void setupOffsetsParameter(Parameter offsetsParameter, List<String> strainNames,
                                       List<Double> strainDates, double earliest, boolean reIndex) {
        offsetsParameter.setDimension(strainNames.size());
        String[] labelArray = new String[strainNames.size()];
        strainNames.toArray(labelArray);
        offsetsParameter.setDimensionNames(labelArray);
        for (int i = 0; i < strainNames.size(); i++) {
            Double offset = strainDates.get(i) - earliest;
            if (offset == null) {
                throw new IllegalArgumentException("Date missing for strain: " + strainNames.get(i));
            }
            offsetsParameter.setParameterValue(reIndex ? tipIndices[i] : i, offset);
        }
        addVariable(offsetsParameter);
    }


    private int[] setupTipIndices(CompoundParameter tipTraitsParameter,
                                 List<String> strainNames) {

        int[] tipIndices = new int[strainNames.size()];
        Arrays.fill(tipIndices, -1);

        for (int i = 0; i < tipTraitsParameter.getParameterCount(); i++) {
            Parameter tip = tipTraitsParameter.getParameter(i);
            String label = tip.getParameterName();
            int index = findStrain(label, strainNames);
            if (index != -1) {
                if (tipIndices[index] != -1) {
                    throw new IllegalArgumentException("Duplicated tip name: " + label);
                }

                tipIndices[index] = i;

                // rather than setting these here, we set them when the locations are set so the changes propagate
                // through to the diffusion model.
            } else {
                // The tree may contain viruses not present in the assay data
                throw new IllegalArgumentException("Unmatched tip name in assay data: " + label + "\nNot yet implemented.");
            }
        }

        return tipIndices;
    }

//    private int[] setupVirusIndices(int[] tipIndices) {
//        int[] virusIndices = new int[tipIndices.length];
//        for (int i = 0; i < tipIndices.length; ++i) {
//            virusIndices[tipIndices[i]] = i;
//        }
//        return virusIndices;
//    }

    private int findStrain(String label, List<String> strainNames) {
        int index = 0;
        for (String strainName : strainNames) {
            if (label.equalsIgnoreCase(strainName)) {
                return index;
            }

            index ++;
        }
        return -1;
    }

    private void setupInitialLocations(double drift) {
        for (int i = 0; i < tipTraitsParameter.getParameterCount(); ++i) {
            double offset = 0.0;
            if (virusOffsetsParameter != null) {
                offset = drift * virusOffsetsParameter.getParameterValue(tipIndices[i]);
            }
            double r = MathUtils.nextGaussian();
            virusSamplingParameter.getParameter(i).setParameterValue(0, r + offset); // TODO Remove
            tipTraitsParameter.getParameter(tipIndices[i]).setParameterValue(tipStartOffset, r + offset);
            if (mdsDimension > 1) {
                for (int j = 1; j < mdsDimension; j++) {
                    r = MathUtils.nextGaussian();
                    virusSamplingParameter.getParameter(i).setParameterValue(j, r); // TODO Remove
                    tipTraitsParameter.getParameter(tipIndices[i]).setParameterValue(tipStartOffset + j, r);
                }
            }
        }
        for (int i = 0; i < serumLocationsParameter.getParameterCount(); i++) {
            double offset = 0.0;
            if (serumOffsetsParameter != null) {
                offset = drift * serumOffsetsParameter.getParameterValue(i);
            }
            double r = MathUtils.nextGaussian() + offset;
            serumLocationsParameter.getParameter(i).setParameterValue(0, r);
            if (mdsDimension > 1) {
                for (int j = 1; j < mdsDimension; j++) {
                    r = MathUtils.nextGaussian();
                    serumLocationsParameter.getParameter(i).setParameterValue(j, r);
                }
            }
        }
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Variable.ChangeType type) {
        if (variable == virusSamplingParameter) { // TODO Remove
            if (index != -1) {
                int loc = index / tipDimension;
//                virusLocationChanged[loc] = true;
                if (tipTraitsParameter != null && tipIndices[loc] != -1) {
                    Parameter location = virusSamplingParameter.getParameter(loc);
                    Parameter tip = tipTraitsParameter.getParameter(tipIndices[loc]);
                    int dim = index % mdsDimension;
                    tip.setParameterValue(dim, location.getParameterValue(dim));
                }
            } else {
//                Arrays.fill(virusLocationChanged, true);
                if (tipTraitsParameter != null) {
                    for (int pIndex = 0; pIndex < virusSamplingParameter.getParameterCount(); ++pIndex) {
                        Parameter location = virusSamplingParameter.getParameter(pIndex);
                        Parameter tip = tipTraitsParameter.getParameter(tipIndices[pIndex]);
                        for (int i = 0; i < tip.getDimension(); ++i) {
                            tip.setParameterValueQuietly(i, location.getParameterValue(i));//
                        }
                    }
                    tipTraitsParameter.fireParameterChangedEvent();
                }

            }

            locationsKnown = false;
        } else if (variable == tipTraitsParameter) {
            locationsKnown = false;
        } else if (variable == serumLocationsParameter) {
            int loc = index / mdsDimension;
//            serumLocationChanged[loc] = true;
            locationsKnown = false;
        } else if (variable == mdsPrecisionParameter) {
//            setLocationChangedFlags(true);
            precisionKnown = false;
        } else if (variable == locationDriftParameter) {
//            setLocationChangedFlags(true);
            locationsKnown = false;
        } else if (variable == virusDriftParameter) {
//            setLocationChangedFlags(true);
            locationsKnown = false;
            throw new IllegalArgumentException("Not yet implemented");
        } else if (variable == serumDriftParameter) {
//            setLocationChangedFlags(true);
            locationsKnown = false;
            throw new IllegalArgumentException("Not yet implemented");
        } else if (variable == serumPotenciesParameter) {
//            serumEffectChanged[index] = true;
            observationsKnown = false;
            throw new IllegalArgumentException("Not yet implemented");
        } else if (variable == serumBreadthsParameter) {
//            serumEffectChanged[index] = true;
            observationsKnown = false;
            throw new IllegalArgumentException("Not yet implemented");
        } else if (variable == virusAviditiesParameter) {
//            virusEffectChanged[index] = true;
            observationsKnown = false;
            throw new IllegalArgumentException("Not yet implemented");
        } else {
            throw new IllegalArgumentException("Not yet implemented");
        }

        fireModelChanged(variable, index);
        likelihoodKnown = false;
    }

    @Override
    protected void storeState() {

        mdsCore.storeState();

        savedLogLikelihood = logLikelihood;
        savedLikelihoodKnown = likelihoodKnown;
        savedLocationsKnown = locationsKnown;
        savedObservationsKnown = observationsKnown;
        savedPrecisionKnown = precisionKnown;
    }

    @Override
    protected void restoreState() {

        mdsCore.restoreState();

        logLikelihood = savedLogLikelihood;
        likelihoodKnown = savedLikelihoodKnown;
        locationsKnown = savedLocationsKnown;
        observationsKnown = savedObservationsKnown;
        precisionKnown = savedPrecisionKnown;
    }

    @Override
    protected void acceptState() {
    }

    public Model getModel() {
        return this;
    }

    public double getLogLikelihood() {
        if (!likelihoodKnown) {
            logLikelihood = computeLogLikelihoodMds();
            likelihoodKnown = true;
        }

        return logLikelihood;
    }

    // This function can be overwritten to implement other sampling densities, i.e. discrete ranks
//    private double computeLogLikelihood() {
//
//        double precision = mdsPrecisionParameter.getParameterValue(0);
//        double sd = 1.0 / Math.sqrt(precision);
//
//        logLikelihood = 0.0;
//        int i = 0;
//
//        double SSE = 0.0;
//        for (Measurement measurement : measurements) {
//
//            if (virusLocationChanged[measurement.virus] || serumLocationChanged[measurement.serum] || virusEffectChanged[measurement.virus] || serumEffectChanged[measurement.serum]) {
//
//                double distance = computeDistance(measurement.virus, measurement.serum);
//                double expectation = calculateBaseline(measurement.virus, measurement.serum) - distance;
//
//                switch (measurement.type) {
//                    case INTERVAL: {
//                        double minTitre = measurement.log2Titre;
//                        double maxTitre = measurement.log2Titre + intervalWidth;
//                        logLikelihoods[i] = computeMeasurementIntervalLikelihood(minTitre, maxTitre, expectation, sd);
//                    } break;
//                    case POINT: {
//                        logLikelihoods[i] = computeMeasurementLikelihood(measurement.log2Titre, expectation, sd);
//                    } break;
//                    case THRESHOLD: {
//                    	if(measurement.isLowerThreshold){
//                    		logLikelihoods[i] = computeMeasurementThresholdLikelihood(measurement.log2Titre, expectation, sd);
//                    	}
//                    	else{
//                    		logLikelihoods[i] = computeMeasurementUpperThresholdLikelihood(measurement.log2Titre, expectation, sd);
//                    	}
//                    } break;
//                    case MISSING:
//                        break;
//                }
//            }
//            logLikelihood += logLikelihoods[i];
//            i++;
//        }
//
//        likelihoodKnown = true;
//
//        setLocationChangedFlags(false);
//        setSerumEffectChangedFlags(false);
//        setVirusEffectChangedFlags(false);
//
//        return logLikelihood;
//    }
//
//    private void setLocationChangedFlags(boolean flag) {
//        for (int i = 0; i < virusLocationChanged.length; i++) {
//            virusLocationChanged[i] = flag;
//        }
//        for (int i = 0; i < serumLocationChanged.length; i++) {
//            serumLocationChanged[i] = flag;
//        }
//    }
//
//    private void setSerumEffectChangedFlags(boolean flag) {
//        for (int i = 0; i < serumEffectChanged.length; i++) {
//            serumEffectChanged[i] = flag;
//        }
//    }
//
//    private void setVirusEffectChangedFlags(boolean flag) {
//        for (int i = 0; i < virusEffectChanged.length; i++) {
//            virusEffectChanged[i] = flag;
//        }
//    }

    // offset virus and serum location when computing
//    protected double computeDistance(int virus, int serum) {
//
//        Parameter vLoc = tipTraitsParameter.getParameter(tipIndices[virus]);
//        Parameter sLoc = serumLocationsParameter.getParameter(serum);
//        double sum = 0.0;
//
//        // first dimension is shifted
//        double vxOffset = 0.0;
//        double sxOffset = 0.0;
//        if (locationDriftParameter != null && virusOffsetsParameter != null && serumOffsetsParameter != null) {
//            vxOffset = locationDriftParameter.getParameterValue(0) * virusOffsetsParameter.getParameterValue(virus);
//            sxOffset = locationDriftParameter.getParameterValue(0) * serumOffsetsParameter.getParameterValue(serum);
//        }
//        if (virusDriftParameter != null && virusOffsetsParameter != null) {
//            vxOffset = virusDriftParameter.getParameterValue(0) * virusOffsetsParameter.getParameterValue(virus);
//        }
//        if (serumDriftParameter != null && serumOffsetsParameter != null) {
//            sxOffset = serumDriftParameter.getParameterValue(0) * serumOffsetsParameter.getParameterValue(serum);
//        }
//
//        double vxLoc = vLoc.getParameterValue(tipStartOffset + 0) + vxOffset;
//        double sxLoc = sLoc.getParameterValue(0) + sxOffset;
//
//        double difference = vxLoc - sxLoc;
//        sum += difference * difference;
//
//        // other dimensions are not
//        for (int i = 1; i < mdsDimension; i++) {
//            difference = vLoc.getParameterValue(tipStartOffset + i) - sLoc.getParameterValue(i);
//            sum += difference * difference;
//        }
//
//        double dist = Math.sqrt(sum);
//
//        if (serumBreadthsParameter != null) {
//            double serumBreadth = serumBreadthsParameter.getParameterValue(serum);
//            dist /= serumBreadth;
//        }
//
//        return dist;
//    }


    // Calculates the expected log2 titre when mapDistance = 0
    private double calculateBaseline(int virus, int serum) {
        double baseline = serumPotenciesParameter.getParameterValue(serum);
        if (virusAviditiesParameter != null) {
            baseline += virusAviditiesParameter.getParameterValue(virus);
        }
        return baseline;
    }

//    private static double computeMeasurementLikelihood(double titre, double expectation, double sd) {
//
//        double lnL = NormalDistribution.logPdf(titre, expectation, sd);
//
//        if (CHECK_INFINITE && Double.isNaN(lnL) || Double.isInfinite(lnL)) {
//            throw new RuntimeException("infinite point measurement");
//        }
//        return lnL;
//    }

//    private static double computeMeasurementThresholdLikelihood(double titre, double expectation, double sd) {
//
//        // real titre is somewhere between -infinity and measured 'titre'
//        // want the lower tail of the normal CDF
//
//        double lnL = NormalDistribution.cdf(titre, expectation, sd, true);          // returns logged CDF
//
//        if (CHECK_INFINITE && Double.isNaN(lnL) || Double.isInfinite(lnL)) {
//            throw new RuntimeException("infinite threshold measurement");
//        }
//        return lnL;
//    }

//    private static double computeMeasurementUpperThresholdLikelihood(double titre, double expectation, double sd) {
//
//        // real titre is somewhere between -infinity and measured 'titre'
//        // want the lower tail of the normal CDF
//    	double L = NormalDistribution.cdf(titre, expectation, sd, false);          // returns  CDF
//    	double lnL = Math.log(1-L);  //get the upper tail probability, then log it
//
//        if (CHECK_INFINITE && Double.isNaN(lnL) || Double.isInfinite(lnL)) {
//            throw new RuntimeException("infinite threshold measurement");
//        }
//        return lnL;
//    }
    
    
//    private static double computeMeasurementIntervalLikelihood(double minTitre, double maxTitre, double expectation, double sd) {
//
//        // real titre is somewhere between measured minTitre and maxTitre
//
//        double cdf1 = NormalDistribution.cdf(maxTitre, expectation, sd, true);     // returns logged CDF
//        double cdf2 = NormalDistribution.cdf(minTitre, expectation, sd, true);     // returns logged CDF
//        double lnL = LogTricks.logDiff(cdf1, cdf2);
//
//        if (CHECK_INFINITE && Double.isNaN(lnL) || Double.isInfinite(lnL)) {
//            // this occurs when the interval is in the far tail of the distribution, cdf1 == cdf2
//            // instead return logPDF of the point
//            lnL = NormalDistribution.logPdf(minTitre, expectation, sd);
//            if (CHECK_INFINITE && Double.isNaN(lnL) || Double.isInfinite(lnL)) {
//                throw new RuntimeException("infinite interval measurement");
//            }
//        }
//        return lnL;
//    }

    public void makeDirty() {

        mdsCore.makeDirty();

        locationsKnown = false;
        observationsKnown = false;
        precisionKnown = false;
        internalGradientKnown = false;
        likelihoodKnown = false;
    }

    public AntigenicGradientWrtParameter wrtFactory(Parameter parameter) {
        if (parameter == tipTraitsParameter) {
            return new AntigenicGradientWrtParameter.VirusLocations(numViruses, numSera, mdsDimension,
                    tipTraitsParameter, layout, tipStartOffset, tipDimension);
        } else if (parameter == serumLocationsParameter) {
            return new AntigenicGradientWrtParameter.SerumLocations(numViruses, numSera, mdsDimension,
                    serumLocationsParameter, layout);
        } else if (parameter == locationDriftParameter) {
            return new AntigenicGradientWrtParameter.Drift(numViruses, numSera, mdsDimension,
                    locationDriftParameter, virusOffsetsParameter, serumOffsetsParameter, layout);
        } else {
            throw new IllegalArgumentException("Not yet implemented");
        }
    }

    private static class Measurement {
        private Measurement(final int virus, final int serum, final double virusDate, final double serumDate, final MeasurementType type, final double titre, final boolean isLowerThreshold) {
            this.virus = virus;
            this.serum = serum;
            this.virusDate = virusDate;
            this.serumDate = serumDate;
            this.type = type;
            this.titre = titre;
            this.log2Titre = Math.log(titre) / Math.log(2);
            this.isLowerThreshold = isLowerThreshold;
        }

        final int virus;
        final int serum;
        final double virusDate;
        final double serumDate;
        final MeasurementType type;
        final double titre;
        final double log2Titre;
        final boolean isLowerThreshold;
    }

    private final List<Measurement> measurements = new ArrayList<>();
    private final List<String> virusNames = new ArrayList<>();
    private final List<String> serumNames = new ArrayList<>();
    private final List<Double> virusDates = new ArrayList<>();
    private final List<Double> serumDates = new ArrayList<>();

    private final int mdsDimension;
    private final int tipDimension;
    private final int tipStartOffset;
    private final double intervalWidth;
    private final Parameter mdsPrecisionParameter;
    private final Parameter locationDriftParameter;
    private final Parameter virusDriftParameter;
    private final Parameter serumDriftParameter;

    private final MatrixParameter virusSamplingParameter;
    private final MatrixParameterInterface serumLocationsParameter;

    private final Parameter virusOffsetsParameter;
    private final Parameter serumOffsetsParameter;

    private final CompoundParameter tipTraitsParameter;
    private final int[] tipIndices;
//    private final int[] virusIndices;

    private final Parameter virusAviditiesParameter;
    private final Parameter serumPotenciesParameter;
    private final Parameter serumBreadthsParameter;

    private double logLikelihood = 0.0;
    private boolean likelihoodKnown = false;

//    private final boolean[] virusLocationChanged;
//    private final boolean[] serumLocationChanged;
//    private final boolean[] serumEffectChanged;
//    private final boolean[] virusEffectChanged;
//    private double[] logLikelihoods;
//    private double[] storedLogLikelihoods;

    @Override
    public Citation.Category getCategory() {
        return Citation.Category.TRAIT_MODELS;
    }

    @Override
    public String getDescription() {
        return "Bayesian Antigenic Cartography framework";
    }

    public List<Citation> getCitations() {
        return Collections.singletonList(CommonCitations.BEDFORD_2015_INTEGRATING);
    }
}

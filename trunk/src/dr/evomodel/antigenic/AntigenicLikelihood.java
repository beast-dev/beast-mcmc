package dr.evomodel.antigenic;

import dr.evolution.util.*;
import dr.inference.model.*;
import dr.math.MathUtils;
import dr.math.LogTricks;
import dr.math.distributions.NormalDistribution;
import dr.util.*;
import dr.xml.*;


import java.io.*;
import java.util.*;
import java.util.logging.Logger;

/**
 * @author Andrew Rambaut
 * @author Trevor Bedford
 * @author Marc Suchard
 * @version $Id$
 */
/*
    Both virus locations and serum locations are shifted by the parameter locationDrift.
    A location is increased by locationDrift x offset.
    Offset is set to 0 for the earliest virus and increasing with difference in date from earliest virus.
    This makes the raw virusLocations and serumLocations parameters not directly interpretable.
*/
public class AntigenicLikelihood extends AbstractModelLikelihood implements Citable {
    private static final boolean CHECK_INFINITE = false;
    private static final boolean USE_THRESHOLDS = true;
    private static final boolean USE_INTERVALS = true;

    public final static String ANTIGENIC_LIKELIHOOD = "antigenicLikelihood";

    // column indices in table
    private static final int VIRUS_ISOLATE = 0;
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

    public AntigenicLikelihood(
            int mdsDimension,
            Parameter mdsPrecisionParameter,
            Parameter locationDriftParameter,
            TaxonList strainTaxa,
            MatrixParameter virusLocationsParameter,
            MatrixParameter serumLocationsParameter,
            MatrixParameter locationsParameter,
            CompoundParameter tipTraitsParameter,
            Parameter virusOffsetsParameter,
            Parameter serumOffsetsParameter,
            Parameter columnParameter,
            Parameter rowParameter,
            DataTable<String[]> dataTable,
            boolean mergeSerumIsolates,
            double intervalWidth) {

        super(ANTIGENIC_LIKELIHOOD);

        List<String> strainNames = new ArrayList<String>();
        List<String> virusNames = new ArrayList<String>();
        List<String> serumNames = new ArrayList<String>();
        Map<String, Double> strainDateMap = new HashMap<String, Double>();

        this.intervalWidth = intervalWidth;
        boolean useIntervals = USE_INTERVALS && intervalWidth > 0.0;

        int thresholdCount = 0;


        earliestDate = Double.POSITIVE_INFINITY;
        for (int i = 0; i < dataTable.getRowCount(); i++) {
            String[] values = dataTable.getRow(i);
            int column = columnLabels.indexOf(values[SERUM_ISOLATE]);
            if (column == -1) {
                columnLabels.add(values[SERUM_ISOLATE]);
                column = columnLabels.size() - 1;
            }

            int columnStrain = -1;
            String columnStrainName;
            if (mergeSerumIsolates) {
                columnStrainName = values[SERUM_STRAIN];
            } else {
                columnStrainName = values[SERUM_ISOLATE];
            }

            if (strainTaxa != null) {
                columnStrain = strainTaxa.getTaxonIndex(columnStrainName);

                throw new UnsupportedOperationException("Should extract dates from taxon list...");
            } else {
                columnStrain = strainNames.indexOf(columnStrainName);
                if (columnStrain == -1) {
                    strainNames.add(columnStrainName);
                    double date = Double.parseDouble(values[SERUM_DATE]);
                    strainDateMap.put(columnStrainName, date);
                    columnStrain = strainNames.size() - 1;
                }
                int thisStrain = serumNames.indexOf(columnStrainName);
                if (thisStrain == -1) {
                    serumNames.add(columnStrainName);
                }
            }

            double columnDate = Double.parseDouble(values[SERUM_DATE]);

            if (columnStrain == -1) {
                throw new IllegalArgumentException("Error reading data table: Unrecognized serum strain name, " + values[SERUM_STRAIN] + ", in row " + (i+1));
            }

            int row = rowLabels.indexOf(values[VIRUS_ISOLATE]);
            if (row == -1) {
                rowLabels.add(values[VIRUS_ISOLATE]);
                row = rowLabels.size() - 1;
            }

            int rowStrain = -1;
            String rowStrainName = values[VIRUS_STRAIN];
            if (strainTaxa != null) {
                rowStrain = strainTaxa.getTaxonIndex(rowStrainName);
            } else {
                rowStrain = strainNames.indexOf(rowStrainName);
                if (rowStrain == -1) {
                    strainNames.add(rowStrainName);
                    double date = Double.parseDouble(values[VIRUS_DATE]);
                    strainDateMap.put(rowStrainName, date);
                    rowStrain = strainNames.size() - 1;
                }
                int thisStrain = virusNames.indexOf(rowStrainName);
                if (thisStrain == -1) {
                    virusNames.add(rowStrainName);
                }
            }
            if (rowStrain == -1) {
                throw new IllegalArgumentException("Error reading data table: Unrecognized virus strain name, " + values[VIRUS_STRAIN] + ", in row " + (i+1));
            }

            double rowDate = Double.parseDouble(values[VIRUS_DATE]);

            boolean isThreshold = false;
            double rawTitre = Double.NaN;
            if (values[TITRE].length() > 0) {
                try {
                    rawTitre = Double.parseDouble(values[TITRE]);
                } catch (NumberFormatException nfe) {
                    // check if threshold below
                    if (values[TITRE].contains("<")) {
                        rawTitre = Double.parseDouble(values[TITRE].replace("<",""));
                        isThreshold = true;
                        thresholdCount++;
                    }
                    // check if threshold above
                    if (values[TITRE].contains(">")) {
                        throw new IllegalArgumentException("Error in measurement: unsupported greater than threshold at row " + (i+1));
                    }
                }
            }

            if (columnDate < earliestDate) {
                earliestDate = columnDate;
            }

            if (rowDate < earliestDate) {
                earliestDate = rowDate;
            }

            MeasurementType type = (isThreshold ? MeasurementType.THRESHOLD : (useIntervals ? MeasurementType.INTERVAL : MeasurementType.POINT));
            Measurement measurement = new Measurement(column, columnStrain, columnDate, row, rowStrain, rowDate, type, rawTitre);

            if (USE_THRESHOLDS || !isThreshold) {
                measurements.add(measurement);
            }

        }

        double[] maxColumnTitre = new double[columnLabels.size()];
        double[] maxRowTitre = new double[rowLabels.size()];
        for (Measurement measurement : measurements) {
            double titre = measurement.log2Titre;
            if (Double.isNaN(titre)) {
                titre = measurement.log2Titre;
            }
            if (titre > maxColumnTitre[measurement.column]) {
                maxColumnTitre[measurement.column] = titre;
            }
            if (titre > maxRowTitre[measurement.row]) {
                maxRowTitre[measurement.row] = titre;
            }
        }

        if (strainTaxa != null) {
            this.strains = strainTaxa;

            // fill in the strain name array for local use
            for (int i = 0; i < strains.getTaxonCount(); i++) {
                strainNames.add(strains.getTaxon(i).getId());
            }

        } else {
            Taxa taxa = new Taxa();
            for (String strain : strainNames) {
                taxa.addTaxon(new Taxon(strain));
            }
            this.strains = taxa;
        }

        this.mdsDimension = mdsDimension;

        this.mdsPrecisionParameter = mdsPrecisionParameter;
        addVariable(mdsPrecisionParameter);

        this.locationDriftParameter = locationDriftParameter;
        addVariable(locationDriftParameter);

        this.locationsParameter = locationsParameter;
        setupLocationsParameter(this.locationsParameter, strainNames);

        this.virusLocationsParameter = virusLocationsParameter;
        if (this.virusLocationsParameter != null) {
            setupSubsetLocationsParameter(virusLocationsParameter, locationsParameter, virusNames);
        }

        this.serumLocationsParameter = serumLocationsParameter;
        if (this.serumLocationsParameter != null) {
            setupSubsetLocationsParameter(serumLocationsParameter, locationsParameter, serumNames);
        }

        this.tipTraitsParameter = tipTraitsParameter;
        if (tipTraitsParameter != null) {
            setupTipTraitsParameter(this.tipTraitsParameter, strainNames);
        }

        this.virusOffsetsParameter = virusOffsetsParameter;
        if (virusOffsetsParameter != null) {
            setupOffsetParameter(virusOffsetsParameter, virusNames, strainDateMap, earliestDate);
        }

        this.serumOffsetsParameter = serumOffsetsParameter;
        if (serumOffsetsParameter != null) {
            setupOffsetParameter(serumOffsetsParameter, serumNames, strainDateMap, earliestDate);
        }

        this.columnEffectsParameter = setupColumnEffectsParameter(columnParameter, maxColumnTitre);

        this.rowEffectsParameter = setupRowEffectsParameter(rowParameter, maxRowTitre);

        StringBuilder sb = new StringBuilder();
        sb.append("\tAntigenicLikelihood:\n");
        sb.append("\t\t" + this.strains.getTaxonCount() + " strains\n");
        sb.append("\t\t" + virusNames.size() + " viruses\n");
        sb.append("\t\t" + serumNames.size() + " sera\n");
        sb.append("\t\t" + columnLabels.size() + " unique columns\n");
        sb.append("\t\t" + rowLabels.size() + " unique rows\n");
        sb.append("\t\t" + measurements.size() + " assay measurements\n");
        if (USE_THRESHOLDS) {
            sb.append("\t\t" + thresholdCount + " thresholded measurements\n");
        }
        if (useIntervals) {
            sb.append("\n\t\tAssuming a log 2 measurement interval width of " + intervalWidth + "\n");
        }
        Logger.getLogger("dr.evomodel").info(sb.toString());

        locationChanged = new boolean[this.locationsParameter.getParameterCount()];
        columnEffectChanged = new boolean[maxColumnTitre.length];
        rowEffectChanged = new boolean[maxRowTitre.length];
        logLikelihoods = new double[measurements.size()];
        storedLogLikelihoods = new double[measurements.size()];

        setupInitialLocations();

        makeDirty();
    }

    private Parameter setupRowEffectsParameter(Parameter rowParameter, double[] maxRowTitre) {
        // If no row parameter is given, then we will only use the column effects
        if (rowParameter != null) {
            rowParameter.addBounds(new Parameter.DefaultBounds(Double.MAX_VALUE, 0.0, 1));
            rowParameter.setDimension(rowLabels.size());
            addVariable(rowParameter);
            String[] labelArray = new String[rowLabels.size()];
            rowLabels.toArray(labelArray);
            rowParameter.setDimensionNames(labelArray);
            for (int i = 0; i < maxRowTitre.length; i++) {
                rowParameter.setParameterValueQuietly(i, maxRowTitre[i]);
            }
        }
        return rowParameter;
    }

    private Parameter setupColumnEffectsParameter(Parameter columnParameter, double[] maxColumnTitre) {
        // If no column parameter is given, make one to hold maximum values for scaling titres...
        if (columnParameter == null) {
            columnParameter = new Parameter.Default("columnEffects");
        } else {
            columnParameter.addBounds(new Parameter.DefaultBounds(Double.MAX_VALUE, 0.0, 1));
            addVariable(columnParameter);
        }

        columnParameter.setDimension(columnLabels.size());
        String[] labelArray = new String[columnLabels.size()];
        columnLabels.toArray(labelArray);
        columnParameter.setDimensionNames(labelArray);
        for (int i = 0; i < maxColumnTitre.length; i++) {
            columnParameter.setParameterValueQuietly(i, maxColumnTitre[i]);
        }

        return columnParameter;
    }

    protected void setupLocationsParameter(MatrixParameter locationsParameter, List<String> strains) {
        locationsParameter.setColumnDimension(mdsDimension);
        locationsParameter.setRowDimension(strains.size());
        for (int i = 0; i < strains.size(); i++) {
            locationsParameter.getParameter(i).setId(strains.get(i));
        }
        addVariable(this.locationsParameter);
    }

    protected void setupSubsetLocationsParameter(MatrixParameter subsetLocationsParameter, MatrixParameter locationsParameter, List<String> strains) {
        for (int i = 0; i < locationsParameter.getParameterCount(); i++) {
            Parameter parameter = locationsParameter.getParameter(i);
            if (strains.contains(parameter.getId())) {
                subsetLocationsParameter.addParameter(parameter);
            }
        }
    }

    private void setupOffsetParameter(Parameter datesParameter, List<String> strainNames, Map<String, Double> strainDateMap, double earliest) {
        datesParameter.setDimension(strainNames.size());
        String[] labelArray = new String[strainNames.size()];
        strainNames.toArray(labelArray);
        datesParameter.setDimensionNames(labelArray);
        for (int i = 0; i < strainNames.size(); i++) {
            Double date = strainDateMap.get(strainNames.get(i)) - new Double(earliest);
            if (date == null) {
                throw new IllegalArgumentException("Date missing for strain: " + strainNames.get(i));
            }
            datesParameter.setParameterValue(i, date);
        }
    }


    private void setupTipTraitsParameter(CompoundParameter tipTraitsParameter, List<String> strainNames) {
        tipIndices = new int[strainNames.size()];

        for (int i = 0; i < strainNames.size(); i++) {
            tipIndices[i] = -1;
        }

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
//                Parameter location = locationsParameter.getParameter(index);
//                for (int dim = 0; dim < mdsDimension; dim++) {
//                    tip.setParameterValue(dim, location.getParameterValue(dim));
//                }
            } else {
                // The tree may contain viruses not present in the assay data
         //       throw new IllegalArgumentException("Unmatched tip name in assay data: " + label);
            }
        }
        // we are only setting this parameter not listening to it:
//        addVariable(this.tipTraitsParameter);
    }

    private final int findStrain(String label, List<String> strainNames) {
        int index = 0;
        for (String strainName : strainNames) {
            if (label.startsWith(strainName)) {
                return index;
            }

            index ++;
        }
        return -1;
    }

    private void setupInitialLocations() {
        for (int i = 0; i < locationsParameter.getParameterCount(); i++) {
            for (int j = 0; j < mdsDimension; j++) {
                double r = MathUtils.nextGaussian();
                locationsParameter.getParameter(i).setParameterValue(j, r);
            }
        }
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Variable.ChangeType type) {
        if (variable == locationsParameter) {
            int loc = index / mdsDimension;
            locationChanged[loc] = true;
            if (tipTraitsParameter != null && tipIndices[loc] != -1) {
                Parameter location = locationsParameter.getParameter(loc);
                Parameter tip = tipTraitsParameter.getParameter(tipIndices[loc]);
                int dim = index % mdsDimension;
                tip.setParameterValue(dim, location.getParameterValue(dim));
            }
        } else if (variable == mdsPrecisionParameter) {
            setLocationChangedFlags(true);
        } else if (variable == locationDriftParameter) {
            setLocationChangedFlags(true);
        } else if (variable == columnEffectsParameter) {
            columnEffectChanged[index] = true;
        } else if (variable == rowEffectsParameter) {
            rowEffectChanged[index] = true;
        } else {
            // could be a derived class's parameter
//            throw new IllegalArgumentException("Unknown parameter");
        }
        likelihoodKnown = false;
    }

    @Override
    protected void storeState() {
        System.arraycopy(logLikelihoods, 0, storedLogLikelihoods, 0, logLikelihoods.length);
    }

    @Override
    protected void restoreState() {
        double[] tmp = logLikelihoods;
        logLikelihoods = storedLogLikelihoods;
        storedLogLikelihoods = tmp;

        likelihoodKnown = false;
    }

    @Override
    protected void acceptState() {
    }

    public Model getModel() {
        return this;
    }

    public double getLogLikelihood() {
        if (!likelihoodKnown) {
            logLikelihood = computeLogLikelihood();
        }

        return logLikelihood;
    }

    // This function can be overwritten to implement other sampling densities, i.e. discrete ranks
    private double computeLogLikelihood() {

        double precision = mdsPrecisionParameter.getParameterValue(0);
        double sd = 1.0 / Math.sqrt(precision);

        logLikelihood = 0.0;
        int i = 0;

        for (Measurement measurement : measurements) {

            if (locationChanged[measurement.rowStrain] || locationChanged[measurement.columnStrain] || columnEffectChanged[measurement.column] || rowEffectChanged[measurement.row]) {

                // the row strain is shifted
                double mapDistance = computeDistance(measurement.rowStrain, measurement.columnStrain, measurement.rowDate, measurement.columnDate);
                double expectation = calculateBaseline(measurement.column, measurement.row) - mapDistance;

                switch (measurement.type) {
                    case INTERVAL: {
                        double minTitre = measurement.log2Titre;
                        double maxTitre = measurement.log2Titre + intervalWidth;
                        logLikelihoods[i] = computeMeasurementIntervalLikelihood(minTitre, maxTitre, expectation, sd);
                    } break;
                    case POINT: {
                        logLikelihoods[i] = computeMeasurementLikelihood(measurement.log2Titre, expectation, sd);
                    } break;
                    case THRESHOLD: {
                        logLikelihoods[i] = computeMeasurementThresholdLikelihood(measurement.log2Titre, expectation, sd);
                    } break;
                    case MISSING:
                        break;
                }
            }
            logLikelihood += logLikelihoods[i];
            i++;
        }

        likelihoodKnown = true;

        setLocationChangedFlags(false);
        setColumnEffectChangedFlags(false);
        setRowEffectChangedFlags(false);

        return logLikelihood;
    }

    private void setLocationChangedFlags(boolean flag) {
        for (int i = 0; i < locationChanged.length; i++) {
            locationChanged[i] = flag;
        }
    }

    private void setColumnEffectChangedFlags(boolean flag) {
        for (int i = 0; i < columnEffectChanged.length; i++) {
            columnEffectChanged[i] = flag;
        }
    }

    private void setRowEffectChangedFlags(boolean flag) {
        for (int i = 0; i < rowEffectChanged.length; i++) {
            rowEffectChanged[i] = flag;
        }
    }

    // offset virus and serum location when computing
    protected double computeDistance(int rowStrain, int columnStrain, double rowDate, double columnDate) {
        if (rowStrain == columnStrain) {
            return 0.0;
        }

        Parameter vLoc = locationsParameter.getParameter(rowStrain);
        Parameter sLoc = locationsParameter.getParameter(columnStrain);
        double sum = 0.0;

        // first dimension is shifted
        double vxOffset = locationDriftParameter.getParameterValue(0) * (rowDate - earliestDate);
        double vxLoc = vLoc.getParameterValue(0) + vxOffset;

        double sxOffset = locationDriftParameter.getParameterValue(0) * (columnDate - earliestDate);
        double sxLoc = sLoc.getParameterValue(0) + sxOffset;

        double difference = vxLoc - sxLoc;
        sum += difference * difference;

        // other dimensions are not
        for (int i = 1; i < mdsDimension; i++) {
            difference = vLoc.getParameterValue(i) - sLoc.getParameterValue(i);
            sum += difference * difference;
        }
        return Math.sqrt(sum);
    }

    /**
     * Calculates the expected log2 titre when mapDistance = 0
     * @param column
     * @param row
     * @return
     */
    private double calculateBaseline(int column, int row) {
        double baseline;
        double columnEffect = columnEffectsParameter.getParameterValue(column);
        if (rowEffectsParameter != null) {
            double rowEffect = rowEffectsParameter.getParameterValue(row);
            baseline = 0.5 * (rowEffect + columnEffect);
        } else {
            baseline = columnEffect;
        }
        return baseline;
    }

    private static double computeMeasurementLikelihood(double titre, double expectation, double sd) {

        double lnL = NormalDistribution.logPdf(titre, expectation, sd);

        if (CHECK_INFINITE && Double.isNaN(lnL) || Double.isInfinite(lnL)) {
            throw new RuntimeException("infinite point measurement");
        }
        return lnL;
    }

    private static double computeMeasurementThresholdLikelihood(double titre, double expectation, double sd) {

        // real titre is somewhere between -infinity and measured 'titre'
        // want the lower tail of the normal CDF

        double lnL = NormalDistribution.cdf(titre, expectation, sd, true);          // returns logged CDF

        if (CHECK_INFINITE && Double.isNaN(lnL) || Double.isInfinite(lnL)) {
            throw new RuntimeException("infinite threshold measurement");
        }
        return lnL;
    }

    private static double computeMeasurementIntervalLikelihood(double minTitre, double maxTitre, double expectation, double sd) {

        // real titre is somewhere between measured minTitre and maxTitre

        double cdf1 = NormalDistribution.cdf(maxTitre, expectation, sd, true);     // returns logged CDF
        double cdf2 = NormalDistribution.cdf(minTitre, expectation, sd, true);     // returns logged CDF
        double lnL = LogTricks.logDiff(cdf1, cdf2);

        if (CHECK_INFINITE && Double.isNaN(lnL) || Double.isInfinite(lnL)) {
            // this occurs when the interval is in the far tail of the distribution, cdf1 == cdf2
            // instead return logPDF of the point
            lnL = NormalDistribution.logPdf(minTitre, expectation, sd);
            if (CHECK_INFINITE && Double.isNaN(lnL) || Double.isInfinite(lnL)) {
                 throw new RuntimeException("infinite interval measurement");
             }
        }
        return lnL;
    }

    public void makeDirty() {
        likelihoodKnown = false;
        setLocationChangedFlags(true);
    }

    private class Measurement {
        private Measurement(final int column, final int columnStrain, final double columnDate, final int row, final int rowStrain, final double rowDate, final MeasurementType type, final double titre) {
            this.column = column;
            this.columnStrain = columnStrain;
            this.columnDate = columnDate;
            this.row = row;
            this.rowStrain = rowStrain;
            this.rowDate = rowDate;

            this.type = type;
            this.titre = titre;
            this.log2Titre = Math.log(titre) / Math.log(2);
        }

        final int column;
        final int row;
        final int columnStrain;
        final int rowStrain;
        final double columnDate;
        final double rowDate;

        final MeasurementType type;
        final double titre;
        final double log2Titre;

    };

    private final List<Measurement> measurements = new ArrayList<Measurement>();
    private final List<String> columnLabels = new ArrayList<String>();
    private final List<String> rowLabels = new ArrayList<String>();


    private final int mdsDimension;
    private final double intervalWidth;
    private final Parameter mdsPrecisionParameter;
    private final Parameter locationDriftParameter;

    private final MatrixParameter locationsParameter;
    private final MatrixParameter virusLocationsParameter;
    private final MatrixParameter serumLocationsParameter;

    private double earliestDate;
    private final Parameter virusOffsetsParameter;
    private final Parameter serumOffsetsParameter;

    private final CompoundParameter tipTraitsParameter;
    private final TaxonList strains;

    private int[] tipIndices;

    private final Parameter columnEffectsParameter;
    private final Parameter rowEffectsParameter;

    private double logLikelihood = 0.0;
    private boolean likelihoodKnown = false;

    private final boolean[] locationChanged;
    private final boolean[] columnEffectChanged;
    private final boolean[] rowEffectChanged;
    private double[] logLikelihoods;
    private double[] storedLogLikelihoods;

// **************************************************************
// XMLObjectParser
// **************************************************************

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {
        public final static String FILE_NAME = "fileName";

        public final static String TIP_TRAIT = "tipTrait";
        public final static String LOCATIONS = "locations";
        public final static String VIRUS_LOCATIONS = "virusLocations";
        public final static String SERUM_LOCATIONS = "serumLocations";
        public static final String MDS_DIMENSION = "mdsDimension";
        public static final String MERGE_SERUM_ISOLATES = "mergeSerumIsolates";
        public static final String INTERVAL_WIDTH = "intervalWidth";
        public static final String MDS_PRECISION = "mdsPrecision";
        public static final String LOCATION_DRIFT = "locationDrift";
        public static final String VIRUS_EFFECTS = "virusEffects";
        public static final String SERUM_EFFECTS = "serumEffects";
        public final static String VIRUS_OFFSETS = "virusOffsets";
        public final static String SERUM_OFFSETS = "serumOffsets";

        public static final String STRAINS = "strains";

        public String getParserName() {
            return ANTIGENIC_LIKELIHOOD;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            String fileName = xo.getStringAttribute(FILE_NAME);
            DataTable<String[]> assayTable;
            try {
                assayTable = DataTable.Text.parse(new FileReader(fileName), true, false);
            } catch (IOException e) {
                throw new XMLParseException("Unable to read assay data from file: " + e.getMessage());
            }
            System.out.println("Loaded HI table file: " + fileName);

            boolean mergeSerumIsolates = xo.getAttribute(MERGE_SERUM_ISOLATES, false);

            int mdsDimension = xo.getIntegerAttribute(MDS_DIMENSION);
            double intervalWidth = 0.0;
            if (xo.hasAttribute(INTERVAL_WIDTH)) {
                intervalWidth = xo.getDoubleAttribute(INTERVAL_WIDTH);
            }

            CompoundParameter tipTraitParameter = null;
            if (xo.hasChildNamed(TIP_TRAIT)) {
                tipTraitParameter = (CompoundParameter) xo.getElementFirstChild(TIP_TRAIT);
            }

            TaxonList strains = null;
            if (xo.hasChildNamed(STRAINS)) {
                strains = (TaxonList) xo.getElementFirstChild(STRAINS);
            }

            MatrixParameter virusLocationsParameter = null;
            if (xo.hasChildNamed(VIRUS_LOCATIONS)) {
                virusLocationsParameter = (MatrixParameter) xo.getElementFirstChild(VIRUS_LOCATIONS);
            }

            MatrixParameter serumLocationsParameter = null;
            if (xo.hasChildNamed(SERUM_LOCATIONS)) {
                serumLocationsParameter = (MatrixParameter) xo.getElementFirstChild(SERUM_LOCATIONS);
            }

            MatrixParameter locationsParameter = (MatrixParameter) xo.getElementFirstChild(LOCATIONS);

            Parameter mdsPrecision = (Parameter) xo.getElementFirstChild(MDS_PRECISION);

            Parameter locationDrift = (Parameter) xo.getElementFirstChild(LOCATION_DRIFT);

            Parameter virusOffsetsParameter = null;
            if (xo.hasChildNamed(VIRUS_OFFSETS)) {
                virusOffsetsParameter = (Parameter) xo.getElementFirstChild(VIRUS_OFFSETS);
            }

            Parameter serumOffsetsParameter = null;
            if (xo.hasChildNamed(SERUM_OFFSETS)) {
                serumOffsetsParameter = (Parameter) xo.getElementFirstChild(SERUM_OFFSETS);
            }

            Parameter columnEffectsParameter = null;
            if (xo.hasChildNamed(SERUM_EFFECTS)) {
                columnEffectsParameter = (Parameter) xo.getElementFirstChild(SERUM_EFFECTS);
            }
            Parameter rowEffectsParameter = null;
            if (xo.hasChildNamed(VIRUS_EFFECTS)) {
                rowEffectsParameter = (Parameter) xo.getElementFirstChild(VIRUS_EFFECTS);
            }

            AntigenicLikelihood AGL = new AntigenicLikelihood(
                    mdsDimension,
                    mdsPrecision,
                    locationDrift,
                    strains,
                    virusLocationsParameter,
                    serumLocationsParameter,
                    locationsParameter,
                    tipTraitParameter,
                    virusOffsetsParameter,
                    serumOffsetsParameter,
                    columnEffectsParameter,
                    rowEffectsParameter,
                    assayTable,
                    mergeSerumIsolates,
                    intervalWidth);

            Logger.getLogger("dr.evomodel").info("Using EvolutionaryCartography model. Please cite:\n" + Utils.getCitationString(AGL));

            return AGL;
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
                AttributeRule.newBooleanRule(MERGE_SERUM_ISOLATES, true, "Should serum isolates that map to the same strain have their locations merged? (defaults to false)"),
                AttributeRule.newDoubleRule(INTERVAL_WIDTH, true, "The width of the titre interval in log 2 space"),
                new ElementRule(STRAINS, TaxonList.class, "A taxon list of strains", true),
                new ElementRule(TIP_TRAIT, CompoundParameter.class, "The parameter of tip locations from the tree", true),
                new ElementRule(LOCATIONS, MatrixParameter.class, "The parameter for locations of all virus and sera"),
                new ElementRule(VIRUS_LOCATIONS, MatrixParameter.class, "The parameter of locations of all virus", true),
                new ElementRule(SERUM_LOCATIONS, MatrixParameter.class, "The parameter of locations of all sera", true),
                new ElementRule(LOCATIONS, MatrixParameter.class),
                new ElementRule(VIRUS_OFFSETS, Parameter.class, "An optional parameter for virus dates to be stored", true),
                new ElementRule(SERUM_OFFSETS, Parameter.class, "An optional parameter for serum dates to be stored", true),
                new ElementRule(SERUM_EFFECTS, Parameter.class, "An optional parameter for column effects", true),
                new ElementRule(VIRUS_EFFECTS, Parameter.class, "An optional parameter for row effects", true),
                new ElementRule(MDS_PRECISION, Parameter.class),
                new ElementRule(LOCATION_DRIFT, Parameter.class)
        };

        public Class getReturnType() {
            return AntigenicLikelihood.class;
        }
    };

    public List<Citation> getCitations() {
        List<Citation> citations = new ArrayList<Citation>();
        citations.add(new Citation(
                new Author[]{
                        new Author("T", "Bedford"),
                        new Author("MA", "Suchard"),
                        new Author("P", "Lemey"),
                        new Author("G", "Dudas"),
                        new Author("V", "Gregory"),
                        new Author("AJ", "Hay"),
                        new Author("JW", "McCauley"),
                        new Author("CA", "Russell"),
                        new Author("DJ", "Smith"),
                        new Author("A", "Rambaut")
                },
                Citation.Status.IN_PREPARATION
        ));
        return citations;
    }

    public static void main(String[] args) {
        double[] titres = {0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0};

        System.out.println("titre\tpoint\tinterval(tail)\tinterval(cdf)\tthreshold");
        for (double titre : titres) {
            double point = AntigenicLikelihood.computeMeasurementLikelihood(titre, 0.0, 1.0);
            double interval = AntigenicLikelihood.computeMeasurementIntervalLikelihood(titre + 1.0, titre, 0.0, 1.0);
            double threshold = AntigenicLikelihood.computeMeasurementThresholdLikelihood(titre, 0.0, 1.0);

            System.out.println(titre + "\t" + point + "\t" + interval + "\t" + threshold);
        }
    }
}

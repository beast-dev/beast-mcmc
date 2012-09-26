package dr.evomodel.antigenic;

import dr.evolution.util.*;
import dr.inference.model.*;
import dr.math.MathUtils;
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
public class AntigenicLikelihood extends AbstractModelLikelihood implements Citable {
    private static final boolean CHECK_INFINITE = false;
    private static final boolean USE_THRESHOLDS = true;
    private static final boolean USE_INTERVALS = true;

    public final static String ANTIGENIC_LIKELIHOOD = "antigenicLikelihood";

    // column indices in table
    private static final int COLUMN_LABEL = 0;
    private static final int SERUM_STRAIN = 2;
    private static final int ROW_LABEL = 1;
    private static final int VIRUS_STRAIN = 3;
    private static final int SERUM_DATE = 4;
    private static final int VIRUS_DATE = 5;
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
            TaxonList strainTaxa,
            MatrixParameter virusLocationsParameter,
            MatrixParameter serumLocationsParameter,
            MatrixParameter locationsParameter,
            CompoundParameter tipTraitsParameter,
            Parameter virusDatesParameter,
            Parameter columnParameter,
            Parameter rowParameter,
            DataTable<String[]> dataTable,
            boolean mergeColumnStrains,
            double intervalWidth,
            List<String> virusLocationStatisticList) {

        super(ANTIGENIC_LIKELIHOOD);

        List<String> strainNames = new ArrayList<String>();
        List<String> virusNames = new ArrayList<String>();
        List<String> serumNames = new ArrayList<String>();
        Map<String, Double> strainDateMap = new HashMap<String, Double>();

        this.intervalWidth = intervalWidth;
        boolean useIntervals = USE_INTERVALS && intervalWidth > 0.0;

        int thresholdCount = 0;

        for (int i = 0; i < dataTable.getRowCount(); i++) {
            String[] values = dataTable.getRow(i);
            int column = columnLabels.indexOf(values[COLUMN_LABEL]);
            if (column == -1) {
                columnLabels.add(values[0]);
                column = columnLabels.size() - 1;
            }

            int columnStrain = -1;
            String columnStrainName;
            if (mergeColumnStrains) {
                columnStrainName = values[SERUM_STRAIN];
            } else {
                columnStrainName = values[COLUMN_LABEL];
            }

            if (strainTaxa != null) {
                columnStrain = strainTaxa.getTaxonIndex(columnStrainName);

                throw new UnsupportedOperationException("Should extract dates from taxon list...");
            } else {
                columnStrain = strainNames.indexOf(columnStrainName);
                if (columnStrain == -1) {
                    strainNames.add(columnStrainName);
                    Double date = Double.parseDouble(values[SERUM_DATE]);
                    strainDateMap.put(columnStrainName, date);
                    columnStrain = strainNames.size() - 1;
                }
                int thisStrain = serumNames.indexOf(columnStrainName);
                if (thisStrain == -1) {
                    serumNames.add(columnStrainName);
                }
            }

            if (columnStrain == -1) {
                throw new IllegalArgumentException("Error reading data table: Unrecognized serum strain name, " + values[SERUM_STRAIN] + ", in row " + (i+1));
            }

            int row = rowLabels.indexOf(values[ROW_LABEL]);
            if (row == -1) {
                rowLabels.add(values[ROW_LABEL]);
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
                    Double date = Double.parseDouble(values[VIRUS_DATE]);
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

            MeasurementType type = (isThreshold ? MeasurementType.THRESHOLD : (useIntervals ? MeasurementType.INTERVAL : MeasurementType.POINT));
            Measurement measurement = new Measurement(column, columnStrain, row, rowStrain, type, rawTitre);

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

        if (virusDatesParameter != null) {
            // this parameter is not used in this class but is setup to be used in other classes
            setupDatesParameter(virusDatesParameter, virusNames, strainDateMap);
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
        logLikelihoods = new double[measurements.size()];
        storedLogLikelihoods = new double[measurements.size()];

        setupInitialLocations(strainNames, strainDateMap);

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

    private void setupDatesParameter(Parameter datesParameter, List<String> strainNames, Map<String, Double> strainDateMap) {
        datesParameter.setDimension(strainNames.size());
        String[] labelArray = new String[strainNames.size()];
        strainNames.toArray(labelArray);
        datesParameter.setDimensionNames(labelArray);
        for (int i = 0; i < strainNames.size(); i++) {
            Double date = strainDateMap.get(strainNames.get(i));
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

    private void setupInitialLocations(List<String> strainNames, Map<String,Double> strainDateMap) {
        double earliestDate = Double.POSITIVE_INFINITY;
        for (double date : strainDateMap.values()) {
            if (earliestDate > date) {
                earliestDate = date;
            }
        }

        for (int i = 0; i < locationsParameter.getParameterCount(); i++) {
            double date = (double) strainDateMap.get(strainNames.get(i));
            double diff = (date-earliestDate);
            locationsParameter.getParameter(i).setParameterValue(0, diff + MathUtils.nextGaussian());

            for (int j = 1; j < mdsDimension; j++) {
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
        } else if (variable == columnEffectsParameter) {
            setLocationChangedFlags(true);
        } else if (variable == rowEffectsParameter) {
            setLocationChangedFlags(true);
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

            if (locationChanged[measurement.rowStrain] || locationChanged[measurement.columnStrain]) {

                double mapDistance = computeDistance(measurement.rowStrain, measurement.columnStrain);
                double logNormalization = calculateTruncationNormalization(mapDistance, sd);

                switch (measurement.type) {
                    case INTERVAL: {
                        // once transformed the lower titre becomes the higher distance
                        double minHiDistance = transformTitre(measurement.log2Titre + 1.0, measurement.column, measurement.row);
                        double maxHiDistance = transformTitre(measurement.log2Titre, measurement.column, measurement.row);
                        logLikelihoods[i] = computeMeasurementIntervalLikelihood(minHiDistance, maxHiDistance, mapDistance, sd) - logNormalization;
                    } break;
                    case POINT: {
                        double hiDistance = transformTitre(measurement.log2Titre, measurement.column, measurement.row);
                        logLikelihoods[i] = computeMeasurementLikelihood(hiDistance, mapDistance, sd) - logNormalization;
                    } break;
                    case THRESHOLD: {
                        double hiDistance = transformTitre(measurement.log2Titre, measurement.column, measurement.row);
                        logLikelihoods[i] = computeMeasurementThresholdLikelihood(hiDistance, mapDistance, sd) - logNormalization;
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

        return logLikelihood;
    }

    private void setLocationChangedFlags(boolean flag) {
        for (int i = 0; i < locationChanged.length; i++) {
            locationChanged[i] = flag;
        }
    }

    protected double computeDistance(int rowStrain, int columnStrain) {
        if (rowStrain == columnStrain) {
            return 0.0;
        }

        Parameter X = locationsParameter.getParameter(rowStrain);
        Parameter Y = locationsParameter.getParameter(columnStrain);
        double sum = 0.0;
        for (int i = 0; i < mdsDimension; i++) {
            double difference = X.getParameterValue(i) - Y.getParameterValue(i);
            sum += difference * difference;
        }
        return Math.sqrt(sum);
    }

    /**
     * Transforms a titre into log2 space and normalizes it with respect to a unit normal
     * @param titre
     * @param column
     * @param row
     * @return
     */
    private double transformTitre(double titre, int column, int row) {
        double t;
        double columnEffect = columnEffectsParameter.getParameterValue(column);
        if (rowEffectsParameter != null) {
            double rowEffect = rowEffectsParameter.getParameterValue(row);
            t = ((rowEffect + columnEffect) * 0.5) - titre;
        } else {
            t = columnEffect - titre;
        }
        return t;
    }

    private static double computeMeasurementIntervalLikelihood(double minDistance, double maxDistance, double mean, double sd) {

        double cdf1 = NormalDistribution.cdf(minDistance, mean, sd, false);
        double cdf2 = NormalDistribution.cdf(maxDistance, mean, sd, false);

        double lnL = Math.log(cdf2 - cdf1);
        if (cdf1 == cdf2) {
            lnL = Math.log(cdf1);
        }
        if (CHECK_INFINITE && Double.isNaN(lnL) || Double.isInfinite(lnL)) {
            throw new RuntimeException("infinite");
        }
        return lnL;
    }

    private static double computeMeasurementIntervalLikelihood_CDF(double minDistance, double maxDistance, double mean, double sd) {

        double cdf1 = NormalDistribution.cdf(minDistance, mean, sd, false);
        double cdf2 = NormalDistribution.cdf(maxDistance, mean, sd, false);

        double lnL = Math.log(cdf1 - cdf2);
        if (cdf1 == cdf2) {
            lnL = Math.log(cdf1);
        }
        if (CHECK_INFINITE && Double.isNaN(lnL) || Double.isInfinite(lnL)) {
            throw new RuntimeException("infinite");
        }
        return lnL;
    }

    private static double computeMeasurementLikelihood(double distance, double mean, double sd) {
        double lnL = NormalDistribution.logPdf(distance, mean, sd);
        if (CHECK_INFINITE && Double.isNaN(lnL) || Double.isInfinite(lnL)) {
            throw new RuntimeException("infinite");
        }
        return lnL;
    }

//    private static double computeMeasurementLowerBoundLikelihood(double transformedMinTitre) {
//        // a lower bound in non-transformed titre so the bottom tail of the distribution
//        double cdf = NormalDistribution.standardTail(transformedMinTitre, false);
//        double lnL = Math.log(cdf);
//        if (CHECK_INFINITE && Double.isNaN(lnL) || Double.isInfinite(lnL)) {
//            throw new RuntimeException("infinite");
//        }
//        return lnL;
//    }

    private static double computeMeasurementThresholdLikelihood(double distance, double mean, double sd) {
        // a upper bound in non-transformed titre so the upper tail of the distribution

        // using special tail function of NormalDistribution (see main() in NormalDistribution for test)
        double tail = NormalDistribution.tailCDF(distance, mean, sd);
        double lnL = Math.log(tail);
        if (CHECK_INFINITE && Double.isNaN(lnL) || Double.isInfinite(lnL)) {
            throw new RuntimeException("infinite");
        }
        return lnL;
    }

    private static double calculateTruncationNormalization(double distance, double sd) {
        return NormalDistribution.cdf(distance, 0.0, sd, true);
    }

    public void makeDirty() {
        likelihoodKnown = false;
        setLocationChangedFlags(true);
    }

    private class Measurement {
        private Measurement(final int column, final int columnStrain, final int row, final int rowStrain, final MeasurementType type, final double titre) {
            this.column = column;
            this.columnStrain = columnStrain;
            this.row = row;
            this.rowStrain = rowStrain;

            this.type = type;
            this.titre = titre;
            this.log2Titre = Math.log(titre) / Math.log(2);
        }

        final int column;
        final int row;
        final int columnStrain;
        final int rowStrain;

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

    private final MatrixParameter locationsParameter;
    private final MatrixParameter virusLocationsParameter;
    private final MatrixParameter serumLocationsParameter;

    private final CompoundParameter tipTraitsParameter;
    private final TaxonList strains;

    private int[] tipIndices;

    private final Parameter columnEffectsParameter;
    private final Parameter rowEffectsParameter;

    private double logLikelihood = 0.0;
    private boolean likelihoodKnown = false;

    private final boolean[] locationChanged;
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
        public final static String VIRUS_DATES = "virusDates";
        public static final String MDS_DIMENSION = "mdsDimension";
        public static final String MERGE_COLUMNS = "mergeColumns";
        public static final String INTERVAL_WIDTH = "intervalWidth";
        public static final String MDS_PRECISION = "mdsPrecision";
        public static final String COLUMN_EFFECTS = "columnEffects";
        public static final String ROW_EFFECTS = "rowEffects";

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

            boolean mergeColumnStrains = xo.getAttribute(MERGE_COLUMNS, false);

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

            Parameter virusDatesParameter = null;
            if (xo.hasChildNamed(VIRUS_DATES)) {
                virusDatesParameter = (Parameter) xo.getElementFirstChild(VIRUS_DATES);
            }

            Parameter mdsPrecision = (Parameter) xo.getElementFirstChild(MDS_PRECISION);

            Parameter columnEffectsParameter = null;
            if (xo.hasChildNamed(COLUMN_EFFECTS)) {
                columnEffectsParameter = (Parameter) xo.getElementFirstChild(COLUMN_EFFECTS);
            }
            Parameter rowEffectsParameter = null;
            if (xo.hasChildNamed(ROW_EFFECTS)) {
                rowEffectsParameter = (Parameter) xo.getElementFirstChild(ROW_EFFECTS);
            }

            AntigenicLikelihood AGL = new AntigenicLikelihood(
                    mdsDimension,
                    mdsPrecision,
                    strains,
                    virusLocationsParameter,
                    serumLocationsParameter,
                    locationsParameter,
                    tipTraitParameter,
                    virusDatesParameter,
                    columnEffectsParameter,
                    rowEffectsParameter,
                    assayTable,
                    mergeColumnStrains,
                    intervalWidth,
                    null);

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
                AttributeRule.newBooleanRule(MERGE_COLUMNS, true, "Should columns with the same strain have their locations merged? (defaults to false)"),
                AttributeRule.newDoubleRule(INTERVAL_WIDTH, true, "The width of the titre interval in log 2 space"),
                new ElementRule(STRAINS, TaxonList.class, "A taxon list of strains", true),
                new ElementRule(TIP_TRAIT, CompoundParameter.class, "The parameter of tip locations from the tree", true),
                new ElementRule(LOCATIONS, MatrixParameter.class, "The parameter for locations of all virus and sera"),
                new ElementRule(VIRUS_LOCATIONS, MatrixParameter.class, "The parameter of locations of all virus", true),
                new ElementRule(SERUM_LOCATIONS, MatrixParameter.class, "The parameter of locations of all sera", true),
                new ElementRule(LOCATIONS, MatrixParameter.class),
                new ElementRule(VIRUS_DATES, Parameter.class, "An optional parameter for strain dates to be stored", true),
                new ElementRule(COLUMN_EFFECTS, Parameter.class, "An optional parameter for column effects", true),
                new ElementRule(ROW_EFFECTS, Parameter.class, "An optional parameter for row effects", true),
                new ElementRule(MDS_PRECISION, Parameter.class)
        };

        public Class getReturnType() {
            return ContinuousAntigenicTraitLikelihood.class;
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
                        new Author("C", "Russell"),
                        new Author("D", "Smith"),
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
            double interval2 = AntigenicLikelihood.computeMeasurementIntervalLikelihood_CDF(titre + 1.0, titre, 0.0, 1.0);
            double threshold = AntigenicLikelihood.computeMeasurementThresholdLikelihood(titre, 0.0, 1.0);

            System.out.println(titre + "\t" + point + "\t" + interval + "\t" + interval2 + "\t" + threshold);
        }
    }
}

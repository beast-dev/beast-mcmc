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
            Parameter serumPotenciesParameter,
            Parameter serumBreadthsParameter,
            Parameter virusAviditiesParameter,
            DataTable<String[]> dataTable,
            boolean mergeIsolates,
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
            int serumIsolate = serumLabels.indexOf(values[SERUM_ISOLATE]);
            if (serumIsolate == -1) {
                serumLabels.add(values[SERUM_ISOLATE]);
                serumIsolate = serumLabels.size() - 1;
            }

            int serumStrain = -1;
            String serumStrainName;
            if (mergeIsolates) {
                serumStrainName = values[SERUM_STRAIN];
            } else {
                serumStrainName = values[SERUM_ISOLATE];
            }

            if (strainTaxa != null) {
                serumStrain = strainTaxa.getTaxonIndex(serumStrainName);

                throw new UnsupportedOperationException("Should extract dates from taxon list...");
            } else {
                serumStrain = strainNames.indexOf(serumStrainName);
                if (serumStrain == -1) {
                    strainNames.add(serumStrainName);
                    double date = Double.parseDouble(values[SERUM_DATE]);
                    strainDateMap.put(serumStrainName, date);
                    serumStrain = strainNames.size() - 1;
                }
                int thisStrain = serumNames.indexOf(serumStrainName);
                if (thisStrain == -1) {
                    serumNames.add(serumStrainName);
                }
            }

            double serumDate = Double.parseDouble(values[SERUM_DATE]);

            if (serumStrain == -1) {
                throw new IllegalArgumentException("Error reading data table: Unrecognized serum strain name, " + values[SERUM_STRAIN] + ", in row " + (i+1));
            }

            int virusIsolate = virusLabels.indexOf(values[VIRUS_ISOLATE]);
            if (virusIsolate == -1) {
                virusLabels.add(values[VIRUS_ISOLATE]);
                virusIsolate = virusLabels.size() - 1;
            }

            int virusStrain = -1;
            String virusStrainName = values[VIRUS_STRAIN];
            if (strainTaxa != null) {
                virusStrain = strainTaxa.getTaxonIndex(virusStrainName);
            } else {
                virusStrain = strainNames.indexOf(virusStrainName);
                if (virusStrain == -1) {
                    strainNames.add(virusStrainName);
                    double date = Double.parseDouble(values[VIRUS_DATE]);
                    strainDateMap.put(virusStrainName, date);
                    virusStrain = strainNames.size() - 1;
                }
                int thisStrain = virusNames.indexOf(virusStrainName);
                if (thisStrain == -1) {
                    virusNames.add(virusStrainName);
                }
            }
            if (virusStrain == -1) {
                throw new IllegalArgumentException("Error reading data table: Unrecognized virus strain name, " + values[VIRUS_STRAIN] + ", in row " + (i+1));
            }

            double virusDate = Double.parseDouble(values[VIRUS_DATE]);

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

            if (serumDate < earliestDate) {
                earliestDate = serumDate;
            }

            if (virusDate < earliestDate) {
                earliestDate = virusDate;
            }

            MeasurementType type = (isThreshold ? MeasurementType.THRESHOLD : (useIntervals ? MeasurementType.INTERVAL : MeasurementType.POINT));
            Measurement measurement = new Measurement(serumIsolate, serumStrain, serumDate, virusIsolate, virusStrain, virusDate, type, rawTitre);

            if (USE_THRESHOLDS || !isThreshold) {
                measurements.add(measurement);
            }

        }

        double[] maxColumnTitre = new double[serumLabels.size()];
        double[] maxRowTitre = new double[virusLabels.size()];
        for (Measurement measurement : measurements) {
            double titre = measurement.log2Titre;
            if (Double.isNaN(titre)) {
                titre = measurement.log2Titre;
            }
            if (titre > maxColumnTitre[measurement.serumIsolate]) {
                maxColumnTitre[measurement.serumIsolate] = titre;
            }
            if (titre > maxRowTitre[measurement.virusIsolate]) {
                maxRowTitre[measurement.virusIsolate] = titre;
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

        this.serumPotenciesParameter = setupSerumPotencies(serumPotenciesParameter, maxColumnTitre);
        this.serumBreadthsParameter = setupSerumBreadths(serumBreadthsParameter);
        this.virusAviditiesParameter = setupVirusAvidities(virusAviditiesParameter, maxRowTitre);

        StringBuilder sb = new StringBuilder();
        sb.append("\tAntigenicLikelihood:\n");
        sb.append("\t\t" + this.strains.getTaxonCount() + " strains\n");
        sb.append("\t\t" + virusNames.size() + " viruses\n");
        sb.append("\t\t" + serumNames.size() + " sera\n");
        sb.append("\t\t" + virusLabels.size() + " unique viruses\n");
        sb.append("\t\t" + serumLabels.size() + " unique sera\n");
        sb.append("\t\t" + measurements.size() + " assay measurements\n");
        if (USE_THRESHOLDS) {
            sb.append("\t\t" + thresholdCount + " thresholded measurements\n");
        }
        if (useIntervals) {
            sb.append("\n\t\tAssuming a log 2 measurement interval width of " + intervalWidth + "\n");
        }
        Logger.getLogger("dr.evomodel").info(sb.toString());

        locationChanged = new boolean[this.locationsParameter.getParameterCount()];
        serumEffectChanged = new boolean[maxColumnTitre.length];
        virusEffectChanged = new boolean[maxRowTitre.length];
        logLikelihoods = new double[measurements.size()];
        storedLogLikelihoods = new double[measurements.size()];

        setupInitialLocations();

        makeDirty();
    }

    private Parameter setupVirusAvidities(Parameter virusAviditiesParameter, double[] maxRowTitre) {
        // If no row parameter is given, then we will only use the serum effects
        if (virusAviditiesParameter != null) {
            virusAviditiesParameter.addBounds(new Parameter.DefaultBounds(Double.MAX_VALUE, 0.0, 1));
            virusAviditiesParameter.setDimension(virusLabels.size());
            addVariable(virusAviditiesParameter);
            String[] labelArray = new String[virusLabels.size()];
            virusLabels.toArray(labelArray);
            virusAviditiesParameter.setDimensionNames(labelArray);
            for (int i = 0; i < maxRowTitre.length; i++) {
                virusAviditiesParameter.setParameterValueQuietly(i, maxRowTitre[i]);
            }
        }
        return virusAviditiesParameter;
    }

    private Parameter setupSerumPotencies(Parameter serumPotenciesParameter, double[] maxColumnTitre) {
        // If no serum potencies parameter is given, make one to hold maximum values for scaling titres...
        if (serumPotenciesParameter == null) {
            serumPotenciesParameter = new Parameter.Default("serumPotencies");
        } else {
            serumPotenciesParameter.addBounds(new Parameter.DefaultBounds(Double.MAX_VALUE, 0.0, 1));
            addVariable(serumPotenciesParameter);
        }

        serumPotenciesParameter.setDimension(serumLabels.size());
        String[] labelArray = new String[serumLabels.size()];
        serumLabels.toArray(labelArray);
        serumPotenciesParameter.setDimensionNames(labelArray);
        for (int i = 0; i < maxColumnTitre.length; i++) {
            serumPotenciesParameter.setParameterValueQuietly(i, maxColumnTitre[i]);
        }

        return serumPotenciesParameter;
    }

    private Parameter setupSerumBreadths(Parameter serumBreadthsParameter) {
        // If no serum breadths parameter is given, then we will only use the serum potencies
        if (serumBreadthsParameter != null) {
            serumBreadthsParameter.addBounds(new Parameter.DefaultBounds(Double.MAX_VALUE, 0.0, 1));
            serumBreadthsParameter.setDimension(serumLabels.size());
            addVariable(serumBreadthsParameter);
            String[] labelArray = new String[serumLabels.size()];
            serumLabels.toArray(labelArray);
            serumBreadthsParameter.setDimensionNames(labelArray);
            for (int i = 0; i < serumLabels.size(); i++) {
                serumBreadthsParameter.setParameterValueQuietly(i, 1.0);
            }
        }
        return serumBreadthsParameter;
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
        } else if (variable == serumPotenciesParameter) {
            serumEffectChanged[index] = true;
        } else if (variable == virusAviditiesParameter) {
            virusEffectChanged[index] = true;
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

            if (locationChanged[measurement.virusStrain] || locationChanged[measurement.serumStrain] || serumEffectChanged[measurement.serumIsolate] || virusEffectChanged[measurement.virusIsolate]) {

                // the row strain is shifted
                double mapDistance = computeDistance(measurement.virusStrain, measurement.serumStrain, measurement.virusDate, measurement.serumDate);
                double expectation = calculateBaseline(measurement.serumIsolate, measurement.virusIsolate) - mapDistance;

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
        setSerumEffectChangedFlags(false);
        setVirusEffectChangedFlags(false);

        return logLikelihood;
    }

    private void setLocationChangedFlags(boolean flag) {
        for (int i = 0; i < locationChanged.length; i++) {
            locationChanged[i] = flag;
        }
    }

    private void setSerumEffectChangedFlags(boolean flag) {
        for (int i = 0; i < serumEffectChanged.length; i++) {
            serumEffectChanged[i] = flag;
        }
    }

    private void setVirusEffectChangedFlags(boolean flag) {
        for (int i = 0; i < virusEffectChanged.length; i++) {
            virusEffectChanged[i] = flag;
        }
    }

    // offset virus and serum location when computing
    protected double computeDistance(int virusStrain, int serumStrain, double virusDate, double serumDate) {
        if (virusStrain == serumStrain) {
            return 0.0;
        }

        Parameter vLoc = locationsParameter.getParameter(virusStrain);
        Parameter sLoc = locationsParameter.getParameter(serumStrain);
        double sum = 0.0;

        // first dimension is shifted
        double vxOffset = locationDriftParameter.getParameterValue(0) * (virusDate - earliestDate);
        double vxLoc = vLoc.getParameterValue(0) + vxOffset;

        double sxOffset = locationDriftParameter.getParameterValue(0) * (serumDate - earliestDate);
        double sxLoc = sLoc.getParameterValue(0) + sxOffset;

        double difference = vxLoc - sxLoc;
        sum += difference * difference;

        // other dimensions are not
        for (int i = 1; i < mdsDimension; i++) {
            difference = vLoc.getParameterValue(i) - sLoc.getParameterValue(i);
            sum += difference * difference;
        }

        double dist = Math.sqrt(sum);
 //       if (serumBreadthsParameter != null) {
 //           double serumBreadth = serumBreadthsParameter.getParameterValue(serumStrain);
 //           dist /= serumBreadth;
 //       }

        return dist;
    }

    /**
     * Calculates the expected log2 titre when mapDistance = 0
     * @param serum
     * @param virus
     * @return
     */
    private double calculateBaseline(int serum, int virus) {
        double baseline;
        double serumEffect = serumPotenciesParameter.getParameterValue(serum);
        if (virusAviditiesParameter != null) {
            double virusEffect = virusAviditiesParameter.getParameterValue(virus);
            baseline = 0.5 * (virusEffect + serumEffect);
        } else {
            baseline = serumEffect;
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
        private Measurement(final int serumIsolate, final int serumStrain, final double serumDate, final int virusIsolate, final int virusStrain, final double virusDate, final MeasurementType type, final double titre) {
            this.serumIsolate = serumIsolate;
            this.serumStrain = serumStrain;
            this.serumDate = serumDate;
            this.virusIsolate = virusIsolate;
            this.virusStrain = virusStrain;
            this.virusDate = virusDate;
            this.type = type;
            this.titre = titre;
            this.log2Titre = Math.log(titre) / Math.log(2);
        }

        final int serumIsolate;
        final int serumStrain;
        final double serumDate;
        final int virusIsolate;
        final int virusStrain;
        final double virusDate;
        final MeasurementType type;
        final double titre;
        final double log2Titre;

    };

    private final List<Measurement> measurements = new ArrayList<Measurement>();
    private final List<String> serumLabels = new ArrayList<String>();
    private final List<String> virusLabels = new ArrayList<String>();


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

    private final Parameter virusAviditiesParameter;
    private final Parameter serumPotenciesParameter;
    private final Parameter serumBreadthsParameter;

    private double logLikelihood = 0.0;
    private boolean likelihoodKnown = false;

    private final boolean[] locationChanged;
    private final boolean[] serumEffectChanged;
    private final boolean[] virusEffectChanged;
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
        public static final String MERGE_ISOLATES = "mergeIsolates";
        public static final String INTERVAL_WIDTH = "intervalWidth";
        public static final String MDS_PRECISION = "mdsPrecision";
        public static final String LOCATION_DRIFT = "locationDrift";
        public static final String VIRUS_AVIDITIES = "virusAvidities";
        public static final String SERUM_POTENCIES = "serumPotencies";
        public static final String SERUM_BREADTHS = "serumBreadths";
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

            boolean mergeIsolates = xo.getAttribute(MERGE_ISOLATES, false);

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

            Parameter serumPotenciesParameter = null;
            if (xo.hasChildNamed(SERUM_POTENCIES)) {
                serumPotenciesParameter = (Parameter) xo.getElementFirstChild(SERUM_POTENCIES);
            }

            Parameter serumBreadthsParameter = null;
            if (xo.hasChildNamed(SERUM_BREADTHS)) {
                serumBreadthsParameter = (Parameter) xo.getElementFirstChild(SERUM_BREADTHS);
            }

            Parameter virusAviditiesParameter = null;
            if (xo.hasChildNamed(VIRUS_AVIDITIES)) {
                virusAviditiesParameter = (Parameter) xo.getElementFirstChild(VIRUS_AVIDITIES);
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
                    serumPotenciesParameter,
                    serumBreadthsParameter,
                    virusAviditiesParameter,
                    assayTable,
                    mergeIsolates,
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
                AttributeRule.newBooleanRule(MERGE_ISOLATES, true, "Should serum isolates that map to the same virus strain have their locations merged? (defaults to false)"),
                AttributeRule.newDoubleRule(INTERVAL_WIDTH, true, "The width of the titre interval in log 2 space"),
                new ElementRule(STRAINS, TaxonList.class, "A taxon list of strains", true),
                new ElementRule(TIP_TRAIT, CompoundParameter.class, "The parameter of tip locations from the tree", true),
                new ElementRule(LOCATIONS, MatrixParameter.class, "The parameter for locations of all virus and sera"),
                new ElementRule(VIRUS_LOCATIONS, MatrixParameter.class, "The parameter of locations of all virus", true),
                new ElementRule(SERUM_LOCATIONS, MatrixParameter.class, "The parameter of locations of all sera", true),
                new ElementRule(LOCATIONS, MatrixParameter.class),
                new ElementRule(VIRUS_OFFSETS, Parameter.class, "An optional parameter for virus dates to be stored", true),
                new ElementRule(SERUM_OFFSETS, Parameter.class, "An optional parameter for serum dates to be stored", true),
                new ElementRule(SERUM_POTENCIES, Parameter.class, "An optional parameter for serum potencies", true),
                new ElementRule(SERUM_BREADTHS, Parameter.class, "An optional parameter for serum breadths", true),
                new ElementRule(VIRUS_AVIDITIES, Parameter.class, "An optional parameter for virus avidities", true),
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

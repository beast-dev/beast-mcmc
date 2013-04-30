package dr.evomodel.antigenic;

import dr.inference.model.*;
import dr.math.LogTricks;
import dr.math.MathUtils;
import dr.math.distributions.NormalDistribution;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;
import dr.util.DataTable;
import dr.xml.*;

import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

/**
 * @author Andrew Rambaut
 * @author Trevor Bedford
 * @author Marc Suchard
 * @version $Id$
 */
/*
    Both virus locations and serum locations are shifted by the parameter locationDrift. A location is increased by
    locationDrift x offset. Offset is set to 0 for the earliest virus and increasing with difference in date from
    earliest virus. This makes the raw virusLocations and serumLocations parameters not directly interpretable.
    Using virus and serum offsets to shift locations.  Offsets can be date of isolation, or can be more general.
*/
/*
    Each HI assay has attached to it, the serum isolate, virus isolate, serum strain and virus strain.  We can
    construct locations by collapsing isolates or by collapsing strains.  Each collapsed point has a serum name and a
    virus name, corresponding to either isolate or strain as desired.  Default to collapsing viruses but not
    collapsing sera.
*/
 /*
    All titres stored in measurements, etc... are converted to log2 while being loaded.
*/
/*
    Ordering is always viruses then sera.
*/
/*
    This is currently not working.  Giving "SEVERE: State was not correctly calculated after an operator move." errors.

 */
public class AntigenicLikelihoodSimple extends AbstractModelLikelihood implements Citable {
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

    public AntigenicLikelihoodSimple(
            int mdsDimension,
            Parameter mdsPrecisionParameter,
            Parameter locationDriftParameter,
            MatrixParameter virusLocationsParameter,
            MatrixParameter serumLocationsParameter,
            MatrixParameter locationsParameter,
            CompoundParameter tipTraitsParameter,
            Parameter virusOffsetsParameter,
            Parameter serumOffsetsParameter,
            Parameter virusEffectsParameter,
            Parameter serumEffectsParameter,
            DataTable<String[]> dataTable,
            boolean mergeVirusIsolates,
            boolean mergeSerumIsolates,
            double intervalWidth) {

        super(ANTIGENIC_LIKELIHOOD);

        List<String> virusNames = new ArrayList<String>();
        List<String> serumNames = new ArrayList<String>();
        List<Double> virusDates = new ArrayList<Double>();
        List<Double> serumDates = new ArrayList<Double>();
        List<Double> virusMaxTitres = new ArrayList<Double>();
        List<Double> serumMaxTitres = new ArrayList<Double>();

        this.intervalWidth = intervalWidth;
        boolean useIntervals = USE_INTERVALS && intervalWidth > 0.0;

        int thresholdCount = 0;
        earliestDate = Double.POSITIVE_INFINITY;

        for (int i = 0; i < dataTable.getRowCount(); i++) {

            String[] values = dataTable.getRow(i);

            // select virus strain or virus isolate as appropriate and add to virus names
            String virusName = "";
            if (mergeVirusIsolates) {
                virusName = values[VIRUS_STRAIN];
            } else {
                virusName = values[VIRUS_ISOLATE];
            }

            int virus = virusNames.indexOf(virusName);
            boolean isNewVirus = false;
            if (virus == -1) {
                isNewVirus = true;
                virusNames.add(virusName);
                virus = virusNames.size() - 1;
            }

            // select serum strain or serum isolate as appropriate and add to serum names
            String serumName = "";
            if (mergeSerumIsolates) {
                serumName = values[SERUM_STRAIN];
            } else {
                serumName = values[SERUM_ISOLATE];
            }

            int serum = serumNames.indexOf(serumName);
            boolean isNewSerum = false;
            if (serum == -1) {
                isNewSerum = true;
                serumNames.add(serumName);
                serum = serumNames.size() - 1;
            }

            // add virus and serum dates
            Double virusDate = new Double(values[VIRUS_DATE]);
            if (isNewVirus) {
                virusDates.add(virusDate);
            }

            if (virusDate.doubleValue() < earliestDate) {
                earliestDate = virusDate.doubleValue();
            }

            Double serumDate = new Double(values[SERUM_DATE]);
            if (isNewSerum) {
                serumDates.add(serumDate);
            }

            if (serumDate.doubleValue() < earliestDate) {
                earliestDate = serumDate.doubleValue();
            }

            // parse and label titre
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
            double titre = Math.log(rawTitre) / Math.log(2);

            // update virus and serum max titres
            if (isNewVirus) {
                Double titreD = new Double(titre);
                virusMaxTitres.add(titreD);
            }
            else {
                double maxTitre = virusMaxTitres.get(virus).doubleValue();
                if (maxTitre < titre) {
                    Double titreD = new Double(titre);
                    virusMaxTitres.set(virus, titreD);
                }
            }

            if (isNewSerum) {
                Double titreD = new Double(titre);
                serumMaxTitres.add(titreD);
            }
            else {
                double maxTitre = serumMaxTitres.get(serum).doubleValue();
                if (maxTitre < titre) {
                    Double titreD = new Double(titre);
                    serumMaxTitres.set(serum, titreD);
                }
            }

            // add measurement
            MeasurementType type = (isThreshold ? MeasurementType.THRESHOLD : (useIntervals ? MeasurementType.INTERVAL : MeasurementType.POINT));
            Measurement measurement = new Measurement(serum, virus, type, titre);
            if (USE_THRESHOLDS || !isThreshold) {
                measurements.add(measurement);
            }

        }

        // setup parameters
        this.mdsDimension = mdsDimension;

        this.mdsPrecisionParameter = mdsPrecisionParameter;
        addVariable(this.mdsPrecisionParameter);

        this.locationDriftParameter = locationDriftParameter;
        addVariable(this.locationDriftParameter);

        this.virusLocationsParameter = virusLocationsParameter;
        setupLocationsParameter(this.virusLocationsParameter, virusNames);
        addVariable(this.virusLocationsParameter);

        this.serumLocationsParameter = serumLocationsParameter;
        setupLocationsParameter(this.serumLocationsParameter, serumNames);
        addVariable(this.serumLocationsParameter);

        this.locationsParameter = locationsParameter;
        setupJointLocationsParameter(this.locationsParameter, this.virusLocationsParameter, this.serumLocationsParameter);
        addVariable(this.locationsParameter);

        this.tipTraitsParameter = tipTraitsParameter;

        this.virusOffsetsParameter = virusOffsetsParameter;
        if (virusOffsetsParameter != null) {
            setupOffsetsParameter(virusOffsetsParameter, virusNames, virusDates, earliestDate);
            addVariable(this.virusOffsetsParameter);
        }

        this.serumOffsetsParameter = serumOffsetsParameter;
        if (serumOffsetsParameter != null) {
            setupOffsetsParameter(serumOffsetsParameter, serumNames, serumDates, earliestDate);
            addVariable(this.serumOffsetsParameter);
        }

        this.virusEffectsParameter = virusEffectsParameter;
        if (virusEffectsParameter != null) {
            setupEffectsParameter(this.virusEffectsParameter, virusNames, virusMaxTitres);
            addVariable(this.virusEffectsParameter);
        }

        this.serumEffectsParameter = serumEffectsParameter;
        setupEffectsParameter(this.serumEffectsParameter, serumNames, serumMaxTitres);
        addVariable(this.serumEffectsParameter);

        StringBuilder sb = new StringBuilder();
        sb.append("\tAntigenicLikelihood:\n");
        sb.append("\t\t" + virusNames.size() + " viruses\n");
        sb.append("\t\t" + serumNames.size() + " sera\n");
        sb.append("\t\t" + measurements.size() + " assay measurements\n");
        if (USE_THRESHOLDS) {
            sb.append("\t\t" + thresholdCount + " thresholded measurements\n");
        }
        if (useIntervals) {
            sb.append("\n\t\tAssuming a log 2 measurement interval width of " + intervalWidth + "\n");
        }
        Logger.getLogger("dr.evomodel").info(sb.toString());

        locationChanged = new boolean[this.locationsParameter.getParameterCount()];
        virusEffectChanged = new boolean[virusNames.size()];
        serumEffectChanged = new boolean[serumNames.size()];
        logLikelihoods = new double[measurements.size()];
        storedLogLikelihoods = new double[measurements.size()];

        makeDirty();
    }

    protected void setupLocationsParameter(MatrixParameter locations, List<String> names) {

        locations.setColumnDimension(mdsDimension);
        locations.setRowDimension(names.size());

        for (int i = 0; i < names.size(); i++) {
            Parameter parameter = locations.getParameter(i);
            parameter.setId(names.get(i));
            for (int j = 0; j < mdsDimension; j++) {
                double r = MathUtils.nextGaussian();
                parameter.setParameterValue(j, r);
            }
        }

    }

    protected void setupJointLocationsParameter(MatrixParameter jointLocations, MatrixParameter virusLocations, MatrixParameter serumLocations) {

        for (int i = 0; i < virusLocations.getParameterCount(); i++) {
            Parameter parameter = virusLocations.getParameter(i);
            jointLocations.addParameter(parameter);
        }

        for (int i = 0; i < serumLocations.getParameterCount(); i++) {
             Parameter parameter = serumLocations.getParameter(i);
             jointLocations.addParameter(parameter);
         }

    }

    protected void setupOffsetsParameter(Parameter offsets, List<String> names, List<Double> dates, double earliest) {

        offsets.setDimension(names.size());

        String[] labelArray = new String[names.size()];
        names.toArray(labelArray);
        offsets.setDimensionNames(labelArray);

        for (int i = 0; i < names.size(); i++) {
            Double offset = dates.get(i) - new Double(earliest);
            if (offset == null) {
                throw new IllegalArgumentException("Date missing for strain: " + names.get(i));
            }
            offsets.setParameterValue(i, offset);
        }

    }

    protected void setupEffectsParameter(Parameter effects, List<String> names, List<Double> maxTitres) {

        effects.setDimension(names.size());

        String[] labelArray = new String[names.size()];
        names.toArray(labelArray);
        effects.setDimensionNames(labelArray);

        for (int i = 0; i < names.size(); i++) {
            Double effect = maxTitres.get(i);
            if (effect == null) {
                throw new IllegalArgumentException("Max titre missing for strain: " + names.get(i));
            }
            effects.setParameterValue(i, effect);
        }

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
        } else if (variable == virusEffectsParameter) {
            virusEffectChanged[index] = true;
        } else if (variable == serumEffectsParameter) {
            serumEffectChanged[index] = true;
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

            if (locationChanged[measurement.virus] || locationChanged[measurement.serum] || virusEffectChanged[measurement.virus] || serumEffectChanged[measurement.serum]) {

                double mapDistance = computeDistance(measurement.virus, measurement.serum);
                double baseline = calculateBaseline(measurement.virus, measurement.serum);
                double expectation = baseline - mapDistance;

                switch (measurement.type) {
                    case INTERVAL: {
                        double minTitre = measurement.titre;
                        double maxTitre = measurement.titre + intervalWidth;
                        logLikelihoods[i] = computeMeasurementIntervalLikelihood(minTitre, maxTitre, expectation, sd);
                    } break;
                    case POINT: {
                        logLikelihoods[i] = computeMeasurementLikelihood(measurement.titre, expectation, sd);
                    } break;
                    case THRESHOLD: {
                        logLikelihoods[i] = computeMeasurementThresholdLikelihood(measurement.titre, expectation, sd);
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
        setVirusEffectChangedFlags(false);
        setSerumEffectChangedFlags(false);

        return logLikelihood;
    }

    private void setLocationChangedFlags(boolean flag) {
        for (int i = 0; i < locationChanged.length; i++) {
            locationChanged[i] = flag;
        }
    }

    private void setVirusEffectChangedFlags(boolean flag) {
        for (int i = 0; i < virusEffectChanged.length; i++) {
            virusEffectChanged[i] = flag;
        }
    }

    private void setSerumEffectChangedFlags(boolean flag) {
        for (int i = 0; i < serumEffectChanged.length; i++) {
            serumEffectChanged[i] = flag;
        }
    }

    // offset virus and serum location when computing
    protected double computeDistance(int virus, int serum) {

        Parameter vLocationParameter = virusLocationsParameter.getParameter(virus);
        Parameter sLocationParameter = serumLocationsParameter.getParameter(serum);

        Parameter sLocation = serumLocationsParameter.getParameter(serum);

        double sum = 0.0;

        // first dimension is shifted
        double vxOffset = locationDriftParameter.getParameterValue(0) * virusOffsetsParameter.getParameterValue(virus);
        double vxLoc = vLocationParameter.getParameterValue(0) + vxOffset;
        double sxOffset = locationDriftParameter.getParameterValue(0) * serumOffsetsParameter.getParameterValue(serum);
        double sxLoc = sLocationParameter.getParameterValue(0) + sxOffset;

        double difference = vxLoc - sxLoc;
        sum += difference * difference;

        // other dimensions are not
        for (int i = 1; i < mdsDimension; i++) {
            difference = vLocationParameter.getParameterValue(i) - sLocationParameter.getParameterValue(i);
            sum += difference * difference;
        }
        return Math.sqrt(sum);
    }

     // Calculates the expected log2 titre when mapDistance = 0
    private double calculateBaseline(int virus, int serum) {
        double baseline;
        double serumEffect = serumEffectsParameter.getParameterValue(serum);
        if (virusEffectsParameter != null) {
            double virusEffect = virusEffectsParameter.getParameterValue(virus);
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
        private Measurement(final int serum, final int virus, final MeasurementType type, final double titre) {
            this.serum = serum;
            this.virus = virus;
            this.type = type;
            this.titre = titre;
        }

        final int serum;
        final int virus;
        final MeasurementType type;
        final double titre;

    };

    private final List<Measurement> measurements = new ArrayList<Measurement>();

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

    private int[] tipIndices;

    private final Parameter virusEffectsParameter;
    private final Parameter serumEffectsParameter;

    private double logLikelihood = 0.0;
    private boolean likelihoodKnown = false;

    private final boolean[] locationChanged;
    private final boolean[] virusEffectChanged;
    private final boolean[] serumEffectChanged;
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
        public static final String MERGE_VIRUS_ISOLATES = "mergeVirusIsolates";
        public static final String MERGE_SERUM_ISOLATES = "mergeSerumIsolates";
        public static final String INTERVAL_WIDTH = "intervalWidth";
        public static final String MDS_PRECISION = "mdsPrecision";
        public static final String LOCATION_DRIFT = "locationDrift";
        public static final String VIRUS_EFFECTS = "virusEffects";
        public static final String SERUM_EFFECTS = "serumEffects";
        public final static String VIRUS_OFFSETS = "virusOffsets";
        public final static String SERUM_OFFSETS = "serumOffsets";

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

            boolean mergeVirusIsolates = xo.getAttribute(MERGE_VIRUS_ISOLATES, true);
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

            MatrixParameter virusLocationsParameter = (MatrixParameter) xo.getElementFirstChild(VIRUS_LOCATIONS);

            MatrixParameter serumLocationsParameter = (MatrixParameter) xo.getElementFirstChild(SERUM_LOCATIONS);

            MatrixParameter locationsParameter = (MatrixParameter) xo.getElementFirstChild(LOCATIONS);

            Parameter mdsPrecision = (Parameter) xo.getElementFirstChild(MDS_PRECISION);

            Parameter locationDrift = null;
            if (xo.hasChildNamed(LOCATION_DRIFT)) {
                locationDrift = (Parameter) xo.getElementFirstChild(LOCATION_DRIFT);
            }

            Parameter virusOffsetsParameter = null;
            if (xo.hasChildNamed(VIRUS_OFFSETS)) {
                virusOffsetsParameter = (Parameter) xo.getElementFirstChild(VIRUS_OFFSETS);
            }

            Parameter serumOffsetsParameter = null;
            if (xo.hasChildNamed(SERUM_OFFSETS)) {
                serumOffsetsParameter = (Parameter) xo.getElementFirstChild(SERUM_OFFSETS);
            }

            Parameter virusEffectsParameter = null;
            if (xo.hasChildNamed(VIRUS_EFFECTS)) {
                virusEffectsParameter = (Parameter) xo.getElementFirstChild(VIRUS_EFFECTS);
            }

            Parameter serumEffectsParameter = (Parameter) xo.getElementFirstChild(SERUM_EFFECTS);

            AntigenicLikelihoodSimple AGL = new AntigenicLikelihoodSimple(
                    mdsDimension,
                    mdsPrecision,
                    locationDrift,
                    virusLocationsParameter,
                    serumLocationsParameter,
                    locationsParameter,
                    tipTraitParameter,
                    virusOffsetsParameter,
                    serumOffsetsParameter,
                    virusEffectsParameter,
                    serumEffectsParameter,
                    assayTable,
                    mergeVirusIsolates,
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
                AttributeRule.newBooleanRule(MERGE_VIRUS_ISOLATES, true, "Should multiple virus isolates from the same strain share antigenic locations? (defaults to true)"),
                AttributeRule.newBooleanRule(MERGE_SERUM_ISOLATES, true, "Should multiple serum isolates from the same strain share antigenic locations? (defaults to false)"),
                AttributeRule.newDoubleRule(INTERVAL_WIDTH, true, "The width of the titre interval in log 2 space"),
                new ElementRule(TIP_TRAIT, CompoundParameter.class, "The parameter of tip locations from the tree", true),
                new ElementRule(LOCATIONS, MatrixParameter.class, "The parameter for locations of all virus and sera"),
                new ElementRule(VIRUS_LOCATIONS, MatrixParameter.class, "The parameter of locations of all viruses"),
                new ElementRule(SERUM_LOCATIONS, MatrixParameter.class, "The parameter of locations of all sera"),
                new ElementRule(VIRUS_OFFSETS, Parameter.class, "An optional parameter for virus offsets to be stored", true),
                new ElementRule(SERUM_OFFSETS, Parameter.class, "An optional parameter for serum offsets to be stored", true),
                new ElementRule(VIRUS_EFFECTS, Parameter.class, "An optional parameter for row effects", true),
                new ElementRule(SERUM_EFFECTS, Parameter.class, "The parameter for serum effects"),
                new ElementRule(MDS_PRECISION, Parameter.class, "The parameter for MDS precision"),
                new ElementRule(LOCATION_DRIFT, Parameter.class, "An optional parameter that shifts locations based on isolation dates", true)
        };

        public Class getReturnType() {
            return AntigenicLikelihoodSimple.class;
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
            double point = AntigenicLikelihoodSimple.computeMeasurementLikelihood(titre, 0.0, 1.0);
            double interval = AntigenicLikelihoodSimple.computeMeasurementIntervalLikelihood(titre + 1.0, titre, 0.0, 1.0);
            double threshold = AntigenicLikelihoodSimple.computeMeasurementThresholdLikelihood(titre, 0.0, 1.0);

            System.out.println(titre + "\t" + point + "\t" + interval + "\t" + threshold);
        }
    }
}

package dr.evomodel.antigenic;

import dr.evolution.util.*;
import dr.inference.model.*;
import dr.math.distributions.NormalDistribution;
import dr.util.*;
import dr.xml.*;

import java.io.*;
import java.util.*;
import java.util.logging.Logger;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class AntigenicLikelihood extends AbstractModelLikelihood implements Citable {

    public final static String ANTIGENIC_LIKELIHOOD = "antigenicLikelihood";

    // column indices in table
    private static final int COLUMN_LABEL = 0;
    private static final int SERUM_STRAIN = 2;
    private static final int ROW_LABEL = 1;
    private static final int VIRUS_STRAIN = 3;
    private static final int SERUM_DATE = 4;
    private static final int VIRUS_DATE = 5;
    private static final int RAW_TITRE = 6;
    private static final int MIN_TITRE = 7;
    private static final int MAX_TITRE = 8;

    public enum MeasurementType {
        INTERVAL,
        POINT,
        UPPER_BOUND,
        LOWER_BOUND,
        MISSING
    }

    public AntigenicLikelihood(
            int mdsDimension,
            Parameter mdsPrecisionParameter,
            TaxonList strainTaxa,
            MatrixParameter locationsParameter,
            Parameter columnEffectsParameter,
            Parameter rowEffectsParameter,
            DataTable<String[]> dataTable,
            List<String> virusLocationStatisticList) {

        super(ANTIGENIC_LIKELIHOOD);

        List<String> strainNames = new ArrayList<String>();
        Map<String, Double> strainDateMap = new HashMap<String, Double>();

        for (int i = 0; i < dataTable.getRowCount(); i++) {
            String[] values = dataTable.getRow(i);
            int column = columnLabels.indexOf(values[COLUMN_LABEL]);
            if (column == -1) {
                columnLabels.add(values[0]);
                column = columnLabels.size() - 1;
            }

            int columnStrain = -1;
            if (strainTaxa != null) {
                columnStrain = strainTaxa.getTaxonIndex(values[SERUM_STRAIN]);
            } else {
                columnStrain = strainNames.indexOf(values[SERUM_STRAIN]);
                if (columnStrain == -1) {
                    strainNames.add(values[SERUM_STRAIN]);
                    columnStrain = strainNames.size() - 1;
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
            if (strainTaxa != null) {
                rowStrain = strainTaxa.getTaxonIndex(values[VIRUS_STRAIN]);
            } else {
                rowStrain = strainNames.indexOf(values[VIRUS_STRAIN]);
                if (rowStrain == -1) {
                    strainNames.add(values[VIRUS_STRAIN]);
                    rowStrain = strainNames.size() - 1;
                }
            }
            if (rowStrain == -1) {
                throw new IllegalArgumentException("Error reading data table: Unrecognized virus strain name, " + values[VIRUS_STRAIN] + ", in row " + (i+1));
            }

            // ignoring strain dates for the moment

            double minTitre = Double.NaN;
            if (values[MIN_TITRE].length() > 0) {
                try {
                    minTitre = Double.parseDouble(values[MIN_TITRE]);
                } catch (NumberFormatException nfe) {
                    // do nothing
                }
            }
            double maxTitre = Double.NaN;
            if (values[MAX_TITRE].length() > 0) {
                try {
                    maxTitre = Double.parseDouble(values[MAX_TITRE]);
                } catch (NumberFormatException nfe) {
                    // do nothing
                }
            }

            MeasurementType type = MeasurementType.INTERVAL;
            if (minTitre == maxTitre) {
                type = MeasurementType.POINT;
            }

            if (Double.isNaN(minTitre) || minTitre == 0.0) {
                if (Double.isNaN(maxTitre)) {
                    throw new IllegalArgumentException("Error in measurement: both min and max titre are at bounds in row " + (i+1));
                }

                type = MeasurementType.UPPER_BOUND;
            } else if (Double.isNaN(maxTitre)) {
                type = MeasurementType.LOWER_BOUND;
            }

                Measurement measurement = new Measurement(column, columnStrain, row, rowStrain, type, minTitre, maxTitre);
            measurements.add(measurement);
        }

        if (strainTaxa != null) {
            this.strains = strainTaxa;
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
        setupLocationsParameter(locationsParameter, strains);
        addVariable(locationsParameter);

        columnEffectsParameter.setDimension(columnLabels.size());
        this.columnEffectsParameter = columnEffectsParameter;
        if (columnEffectsParameter != null) {
            addVariable(columnEffectsParameter);
        }

        rowEffectsParameter.setDimension(rowLabels.size());
        this.rowEffectsParameter = rowEffectsParameter;
        if (rowEffectsParameter != null) {
            addVariable(rowEffectsParameter);
        }

        StringBuilder sb = new StringBuilder();
        sb.append("\tAntigenicLikelihood:\n");
        sb.append("\t\t" + this.strains.getTaxonCount() + " strains\n");
        sb.append("\t\t" + columnLabels.size() + " unique columns\n");
        sb.append("\t\t" + rowLabels.size() + " unique rows\n");
        sb.append("\t\t" + measurements.size() + " assay measurements\n");
        Logger.getLogger("dr.evomodel").info(sb.toString());
    }

    protected void setupLocationsParameter(MatrixParameter locationsParameter, TaxonList strains) {
        locationsParameter.setColumnDimension(mdsDimension);
        locationsParameter.setRowDimension(strains.getTaxonCount());
        for (int i = 0; i < strains.getTaxonCount(); i++) {
            locationsParameter.getParameter(i).setId(strains.getTaxon(i).getId());
        }
    }


    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Variable.ChangeType type) {
        if (variable == locationsParameter) {
        } else if (variable == mdsPrecisionParameter) {
        } else if (variable == columnEffectsParameter) {
        } else if (variable == rowEffectsParameter) {
        } else {
            // could be a derived class's parameter
//            throw new IllegalArgumentException("Unknown parameter");
        }
        likelihoodKnown = false;
    }

    @Override
    protected void storeState() {
    }

    @Override
    protected void restoreState() {
        likelihoodKnown = false;
    }

    @Override
    protected void acceptState() {
    }

    @Override
    public Model getModel() {
        return null;
    }

    @Override
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
        for (Measurement measurement : measurements) {
            double distance = computeDistance(measurement.rowStrain, measurement.columnStrain);

            double logNormalization = calculateTruncationNormalization(distance, sd);

            switch (measurement.type) {
                case INTERVAL: {
                    double minTitre = transformTitre(measurement.minTitre, measurement.column, measurement.row);
                    double maxTitre = transformTitre(measurement.maxTitre, measurement.column, measurement.row);
                    logLikelihood += computeMeasurementIntervalLikelihood(distance, minTitre, maxTitre, sd) - logNormalization;
                } break;
                case POINT: {
                    double titre = transformTitre(measurement.minTitre, measurement.column, measurement.row);
                    logLikelihood += computeMeasurementLikelihood(distance, titre, sd) - logNormalization;
                } break;
                case LOWER_BOUND: {
                    double minTitre = transformTitre(measurement.minTitre, measurement.column, measurement.row);
                    logLikelihood += computeMeasurementUpperTailLikelihood(distance, minTitre, sd) - logNormalization;
                } break;
                case UPPER_BOUND: {
                    double maxTitre = transformTitre(measurement.maxTitre, measurement.column, measurement.row);
                    logLikelihood += computeMeasurementLowerTailLikelihood(distance, maxTitre, sd) - logNormalization;
                } break;
                case MISSING:
                    break;
            }
        }

        likelihoodKnown = true;

        return logLikelihood;
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

    private double transformTitre(double titre, int column, int row) {
        double rowEffect = rowEffectsParameter.getParameterValue(row);
        double columnEffect = columnEffectsParameter.getParameterValue(column);

        double t = ((rowEffect + columnEffect) * 0.5) - (Math.log(titre) / Math.log(2));
        return t;
    }

    private double computeMeasurementIntervalLikelihood(double distance, double minTitre, double maxTitre, double sd) {
        double cdf1 = NormalDistribution.cdf(minTitre, distance, sd, false);
        double cdf2 = NormalDistribution.cdf(maxTitre, distance, sd, false);

        double lnL = Math.log(cdf1 - cdf2);
         if (cdf1 == cdf2) {
            lnL = Math.log(cdf1);
        }
        if (Double.isNaN(lnL) || Double.isInfinite(lnL)) {
            throw new RuntimeException("infinite");
        }
        return lnL;
    }

    private double computeMeasurementLikelihood(double distance, double titre, double sd) {
        double lnL = Math.log(NormalDistribution.pdf(titre, distance, sd));
        if (Double.isNaN(lnL) || Double.isInfinite(lnL)) {
            throw new RuntimeException("infinite");
        }
        return lnL;
 }

    private double computeMeasurementLowerTailLikelihood(double distance, double minTransformedTitre, double sd) {
        // using special tail function of NormalDistribution (see main() in NormalDistribution for test)
        double tail = NormalDistribution.tailCDF(minTransformedTitre, distance, sd);
        double lnL = Math.log(tail);
        if (Double.isNaN(lnL) || Double.isInfinite(lnL)) {
            throw new RuntimeException("infinite");
        }
        return lnL;
   }

    private double computeMeasurementUpperTailLikelihood(double distance, double maxTransformedTitre, double sd) {
        double cdf = NormalDistribution.cdf(maxTransformedTitre, distance, sd, false);
        double lnL = Math.log(cdf);
        if (Double.isNaN(lnL) || Double.isInfinite(lnL)) {
            throw new RuntimeException("infinite");
        }
        return lnL;
   }

    private double calculateTruncationNormalization(double distance, double sd) {
        return NormalDistribution.cdf(distance, 0.0, sd, true);
    }

    @Override
    public void makeDirty() {
    }

    private class Measurement {
        private Measurement(final int column, final int columnStrain, final int row, final int rowStrain, final MeasurementType type, final double minTitre, final double maxTitre) {
            this.column = column;
            this.columnStrain = columnStrain;
            this.row = row;
            this.rowStrain = rowStrain;

            this.type = type;
            this.minTitre = minTitre;
            this.maxTitre = maxTitre;
        }

        final int column;
        final int row;
        final int columnStrain;
        final int rowStrain;

        final MeasurementType type;
        final double minTitre;
        final double maxTitre;

    };

    private final List<Measurement> measurements = new ArrayList<Measurement>();
    private final List<String> columnLabels = new ArrayList<String>();
    private final List<String> rowLabels = new ArrayList<String>();


    private final int mdsDimension;
    private final Parameter mdsPrecisionParameter;
    private final MatrixParameter locationsParameter;
    private final TaxonList strains;
    //    private final CompoundParameter tipTraitParameter;
    private final Parameter columnEffectsParameter;
    private final Parameter rowEffectsParameter;

    private double logLikelihood = 0.0;
    private boolean likelihoodKnown = false;

    // **************************************************************
    // XMLObjectParser
    // **************************************************************

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {
        public final static String FILE_NAME = "fileName";

        public final static String TIP_TRAIT = "tipTrait";
        public final static String LOCATIONS = "locations";
        public static final String MDS_DIMENSION = "mdsDimension";
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

            int mdsDimension = xo.getIntegerAttribute(MDS_DIMENSION);

//            CompoundParameter tipTraitParameter = null;
//            if (xo.hasChildNamed(TIP_TRAIT)) {
//                tipTraitParameter = (CompoundParameter) xo.getElementFirstChild(TIP_TRAIT);
//            }

            TaxonList strains = null;
            if (xo.hasChildNamed(STRAINS)) {
                strains = (TaxonList) xo.getElementFirstChild(STRAINS);
            }

            MatrixParameter locationsParameter = (MatrixParameter) xo.getElementFirstChild(LOCATIONS);

            Parameter mdsPrecision = (Parameter) xo.getElementFirstChild(MDS_PRECISION);

            Parameter columnEffectsParameter = (Parameter) xo.getElementFirstChild(COLUMN_EFFECTS);

            Parameter rowEffectsParameter = (Parameter) xo.getElementFirstChild(ROW_EFFECTS);

            AntigenicLikelihood AGL = new AntigenicLikelihood(
                    mdsDimension,
                    mdsPrecision,
                    strains,
                    locationsParameter,
                    columnEffectsParameter,
                    rowEffectsParameter,
                    assayTable,
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
                new ElementRule(STRAINS, TaxonList.class, "A taxon list of strains", true),
//                new ElementRule(TIP_TRAIT, CompoundParameter.class, "The parameter of tip locations from the tree", true),
                new ElementRule(LOCATIONS, MatrixParameter.class),
                new ElementRule(COLUMN_EFFECTS, Parameter.class),
                new ElementRule(ROW_EFFECTS, Parameter.class),
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
}

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
            TaxonList strains,
            CompoundParameter tipTraitParameter,
            MatrixParameter locationsParameter,
            Parameter columnEffectsParameter,
            Parameter rowEffectsParameter,
            DataTable<String[]> dataTable,
            List<String> virusLocationStatisticList) {

        super(ANTIGENIC_LIKELIHOOD);

        this.mdsDimension = mdsDimension;
        this.mdsPrecisionParameter = mdsPrecisionParameter;
        addVariable(mdsPrecisionParameter);

        this.locationsParameter = locationsParameter;
        addVariable(locationsParameter);

        this.strains = strains;

        this.tipTraitParameter = tipTraitParameter;
        if (tipTraitParameter != null) {
            addVariable(tipTraitParameter);
        }

        this.columnEffectsParameter = columnEffectsParameter;
        if (columnEffectsParameter != null) {
            addVariable(columnEffectsParameter);
        }

        this.rowEffectsParameter = rowEffectsParameter;
        if (rowEffectsParameter != null) {
            addVariable(rowEffectsParameter);
        }

        for (int i = 0; i < dataTable.getRowCount(); i++) {
            String[] values = dataTable.getRow(i);
            int column = columnLabels.indexOf(values[0]);
            if (column == -1) {
                columnLabels.add(values[0]);
                column = columnLabels.size() - 1;
            }

            int columnStrain = strains.getTaxonIndex(values[1]);
            if (columnStrain == -1) {
                throw new IllegalArgumentException("Error reading data table: Unrecognized column strain name, " + values[1] + ", in row " + (i+1));
            }

            int row = rowLabels.indexOf(values[2]);
            if (row == -1) {
                rowLabels.add(values[2]);
                row = rowLabels.size() - 1;
            }

            int rowStrain = strains.getTaxonIndex(values[3]);
            if (rowStrain == -1) {
                throw new IllegalArgumentException("Error reading data table: Unrecognized row strain name, " + values[3] + ", in row " + (i+1));
            }

            double minTitre = Double.NaN;
            if (values[4].length() > 0) {
                try {
                    minTitre = Double.parseDouble(values[4]);
                } catch (NumberFormatException nfe) {
                    // do nothing
                }
            }
            double maxTitre = Double.NaN;
            if (values[5].length() > 0) {
                try {
                    maxTitre = Double.parseDouble(values[5]);
                } catch (NumberFormatException nfe) {
                    // do nothing
                }
            }

            MeasurementType type = MeasurementType.INTERVAL;
            if (minTitre == maxTitre) {
                type = MeasurementType.POINT;
            }

            if (Double.isNaN(maxTitre)) {
                type = MeasurementType.UPPER_BOUND;
            } else if (Double.isNaN(minTitre)) {
                type = MeasurementType.LOWER_BOUND;
            }

            Measurement measurement = new Measurement(column, columnStrain, row, rowStrain, type, minTitre, maxTitre);
            measurements.add(measurement);
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
                case UPPER_BOUND: {
                    double minTitre = transformTitre(measurement.minTitre, measurement.column, measurement.row);
                    logLikelihood += computeMeasurementUpperTailLikelihood(distance, minTitre, sd) - logNormalization;
                } break;
                case LOWER_BOUND: {
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

        return (Math.log(Math.sqrt(rowEffect * columnEffect)) - Math.log(titre)) / Math.log(2);
    }

    private double computeMeasurementIntervalLikelihood(double distance, double minTitre, double maxTitre, double sd) {
        double cdf1 = NormalDistribution.cdf(minTitre, distance, sd, false);
        double cdf2 = NormalDistribution.cdf(maxTitre, distance, sd, false);

        return Math.log(cdf2 - cdf1);
    }

    private double computeMeasurementLikelihood(double distance, double titre, double sd) {
        return Math.log(NormalDistribution.pdf(titre, distance, sd));
    }

    private double computeMeasurementLowerTailLikelihood(double distance, double maxTitre, double sd) {
        double cdf = NormalDistribution.cdf(maxTitre, distance, sd, false);
        return Math.log(cdf);
    }

    private double computeMeasurementUpperTailLikelihood(double distance, double minTitre, double sd) {
//        double cdf = NormalDistribution.cdf(maxTitre, distance, sd, false);
//        double tail = 1.0 - cdf;
        // using special tail function of NormalDistribution (see main() in NormalDistribution for test)
        double tail = NormalDistribution.tailCDF(minTitre, distance, sd);
        return Math.log(tail);
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
    private final CompoundParameter tipTraitParameter;
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
                assayTable = DataTable.Text.parse(new FileReader(fileName));
            } catch (IOException e) {
                throw new XMLParseException("Unable to read assay data from file: " + e.getMessage());
            }

            int mdsDimension = xo.getIntegerAttribute(MDS_DIMENSION);

            CompoundParameter tipTraitParameter = null;
            if (xo.hasChildNamed(TIP_TRAIT)) {
                tipTraitParameter = (CompoundParameter) xo.getElementFirstChild(TIP_TRAIT);
            }

            TaxonList strains = (TaxonList) xo.getElementFirstChild(STRAINS);

            MatrixParameter locationsParameter = (MatrixParameter) xo.getElementFirstChild(LOCATIONS);

            Parameter mdsPrecision = (Parameter) xo.getElementFirstChild(MDS_PRECISION);

            Parameter columnEffectsParameter = (Parameter) xo.getElementFirstChild(COLUMN_EFFECTS);

            Parameter rowEffectsParameter = (Parameter) xo.getElementFirstChild(ROW_EFFECTS);

            AntigenicLikelihood AGL = new AntigenicLikelihood(
                    mdsDimension,
                    mdsPrecision,
                    strains,
                    tipTraitParameter,
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
                new ElementRule(TIP_TRAIT, CompoundParameter.class, "The parameter of tip locations from the tree", true),
                new ElementRule(LOCATIONS, MatrixParameter.class),
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
                        new Author("C", "Russell"),
                        new Author("D", "Smith"),
                        new Author("A", "Rambaut")
                },
                Citation.Status.IN_PREPARATION
        ));
        return citations;
    }
}

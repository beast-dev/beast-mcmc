package dr.evomodel.antigenic;

import dr.evolution.util.*;
import dr.inference.model.*;
import dr.util.*;
import dr.xml.*;

import java.io.*;
import java.util.*;
import java.util.logging.Logger;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class AntigenicLikelihood implements Citable {

    public final static String ANTIGENIC_LIKELIHOOD = "antigenicLikelihood";

    public AntigenicLikelihood(
            int mdsDimension,
            Parameter mdsPrecision,
            TaxonList strains,
            CompoundParameter tipTraitParameter,
            MatrixParameter locationsParameter,
            Parameter columnEffectsParameter,
            Parameter rowEffectsParameter,
            DataTable<String[]> dataTable,
            List<String> virusLocationStatisticList) {

        for (int i = 0; i < dataTable.getRowCount(); i++) {
            String[] values = dataTable.getRow(i);
            columnLabels.add(values[0]);
            int column = columnLabels.size() - 1;

            int columnStrain = strains.getTaxonIndex(values[1]);
            if (columnStrain == -1) {
                throw new IllegalArgumentException("Error reading data table: Unrecognized column strain name, " + values[1] + ", in row " + (i+1));
            }

            rowLabels.add(values[2]);
            int row = rowLabels.size() - 1;

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

            Measurement measurement = new Measurement(column, columnStrain, row, rowStrain, minTitre, maxTitre);
            measurements.add(measurement);
        }
    }

    private class Measurement {
        private Measurement(final int column, final int columnStrain, final int row, final int rowStrain, final double minTitre, final double maxTitre) {
            this.column = column;
            this.columnStrain = columnStrain;
            this.minTitre = minTitre;
            this.maxTitre = maxTitre;
            this.row = row;
            this.rowStrain = rowStrain;
        }

        final int column;
        final int row;
        final int columnStrain;
        final int rowStrain;

        final double minTitre;
        final double maxTitre;

    };

    private List<Measurement> measurements = new ArrayList<Measurement>();
    private List<String> columnLabels = new ArrayList<String>();
    private List<String> rowLabels = new ArrayList<String>();

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

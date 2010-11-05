package dr.evomodel.continuous;

import dr.evolution.util.Taxon;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.treelikelihood.TreeTraitParserUtilities;
import dr.inference.model.*;
import dr.math.distributions.NormalDistribution;
import dr.util.DataTable;
import dr.xml.*;

import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

/**
 * @author Andrew Rambaut
 * @author Marc Suchard
 * @version $Id$
 */
public class AntigenicTraitLikelihood extends AbstractModelLikelihood {

    public final static String ANTIGENIC_TRAIT_LIKELIHOOD = "antigenicTraitLikelihood";

    public final static String FILE_NAME = "fileName";

    public final static String TIP_TRAIT = "tipTrait";
    public final static String VIRUS_LOCATIONS = "virusLocations";
    public final static String SERUM_LOCATIONS = "serumLocations";
    public static final String MDS_PRECISION = "mdsPrecision";

    public AntigenicTraitLikelihood(
            TreeModel tree,
            Parameter mdsPrecision,
            CompoundParameter tipTraitParameter,
            CompoundParameter virusLocationsParameter,
            CompoundParameter serumLocationsParameter,
            DataTable<double[]> dataTable) {

        super(ANTIGENIC_TRAIT_LIKELIHOOD);

        String[] virusNames = dataTable.getRowLabels();
        String[] serumNames = dataTable.getColumnLabels();

//        int tipCount = tipTraitParameter.getDimension();

        // the total number of viruses is the number of rows in the table
        int virusCount = dataTable.getRowCount();
        // the number of sera is the number of columns
        int serumCount = dataTable.getColumnCount();

        // the tip -> virus map
        tipIndices = new int[tree.getExternalNodeCount()];

        Map<String, Integer> tipNameMap = new HashMap<String, Integer>();
        for (int i = 0; i < tree.getExternalNodeCount(); i++) {
            Taxon taxon = tree.getNodeTaxon(tree.getExternalNode(i));
            tipNameMap.put(taxon.getId(), i);

            tipIndices[i] = -1;
        }


        // the virus -> tip map
        virusIndices = new int[virusCount];

        // a set of vectors for each virus giving serum indices for which assay data is available
        measuredSerumIndices = new int[virusCount][];

        // a compressed (no missing values) set of measured assay values between virus and sera.
        this.assayTable = new double[virusCount][];

        // a cache of virus to serum distances (serum indices given by array above).
        cachedDistances = new double[virusCount][];

        int totalMeasurementCount = 0;
        for (int i = 0; i < virusCount; i++) {
            virusIndices[i] = -1;

            double[] dataRow = dataTable.getRow(i);

            // if the virus is in the tree then add a entry to map tip to virus
            Integer tipIndex = tipNameMap.get(virusNames[i]);
            if (tipIndex != null) {
                tipIndices[tipIndex] = i;
                virusIndices[i] = tipIndex;
            } else {
                System.err.println("Virus, " + virusNames[i] + ", not found in tree");
            }

            int measuredCount = 0;
            for (int j = 0; j < serumCount; i++) {
                if (!Double.isNaN(dataRow[j]) && dataRow[j] > 0) {
                    measuredCount ++;
                }
            }

            assayTable[i] = new double[measuredCount];
            measuredSerumIndices[i] = new int[measuredCount];
            cachedDistances[i] = new double[measuredCount];

            int k = 0;
            for (int j = 0; j < serumCount; i++) {
                if (!Double.isNaN(dataRow[j]) && dataRow[j] > 0) {
                    this.assayTable[i][k] = transform(dataRow[k]);
                    measuredSerumIndices[i][k] = j;
                    k ++;
                }
            }
            totalMeasurementCount += measuredCount;
        }

        this.totalMeasurementCount = totalMeasurementCount;

        for (int i = 0; i < tree.getExternalNodeCount(); i++) {
            Taxon taxon = tree.getNodeTaxon(tree.getExternalNode(i));
            if (tipIndices[i] == -1) {
                System.err.println("Tree tip, " + taxon.getId() + ", not found in virus assay table");
            }
        }

        this.tipTraitParameter = tipTraitParameter;

        this.virusLocationsParameter = virusLocationsParameter;
        virusLocationsParameter.setDimension(virusCount);
        addVariable(virusLocationsParameter);

        this.serumLocationsParameter = serumLocationsParameter;
        serumLocationsParameter.setDimension(serumCount);
        addVariable(serumLocationsParameter);

        // we don't listen to the tip trait parameter because we are setting that

        this.mdsParameter = mdsPrecision;
        addVariable(mdsPrecision);

        this.isLeftTruncated = true; // Re-normalize likelihood for strictly positive distances
    }

    private double transform(final double value) {
        // transform to log_2
        return Math.log(value) / Math.log(2.0);
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Variable.ChangeType type) {
        // TODO Flag which cachedDistances or mdsPrecision need to be updated
    }

    @Override
    protected void storeState() {
    }

    @Override
    protected void restoreState() {
    }

    @Override
    protected void acceptState() {
        // do nothing
    }

    public Model getModel() {
        return this;
    }

    public double getLogLikelihood() {
        // TODO Only recompute when not known: distances or mdsPrecision changed
        return calculateLogLikelihood();
    }

    private double calculateLogLikelihood() {

        // TODO Only recompute if distances changed
        calculateDistances(); // There is no easy way to loop through non-missing entries of double[][] cachedDistances

        return computeLogLikelihood();
    }

    // This function can be overwritten to implement other sampling densities, i.e. discrete ranks
    protected double computeLogLikelihood() {
        double SSR = calculateSumOfSquaredResiduals();

        double precision = mdsParameter.getParameterValue(0);
        double logLikelihood = (totalMeasurementCount / 2) * Math.log(precision) - 2.0 * precision * SSR;

        if (isLeftTruncated) {
            logLikelihood -= calculateTruncation(precision);
        }
        return logLikelihood;
    }

    private double calculateTruncation(double precision) {
        double sum = 0.0;
        double sd = 1.0 / Math.sqrt(precision);
        for (int i = 0; i < assayTable.length; i++) {
            for (int j = 0; j < assayTable[i].length; j++) {
                sum += Math.log(NormalDistribution.cdf(cachedDistances[i][j], 0.0, sd));
            }
        }
        return sum;
    }

    private double calculateSumOfSquaredResiduals() {
        double sum = 0.0;
        for (int i = 0; i < assayTable.length; i++) {
            for (int j = 0; j < assayTable[i].length; j++) {
                double residual = cachedDistances[i][j] - assayTable[i][j];
                sum += residual * residual;
            }
        }
        return sum;
    }

    private void calculateDistances() {
        for (int i = 0; i < assayTable.length; i++) {
            for (int j = 0; j < assayTable[i].length; j++) {
                cachedDistances[i][j] = calculateDistance(virusLocationsParameter.getParameter(i),
                        serumLocationsParameter.getParameter(measuredSerumIndices[i][j]));
            }
        }
    }

    private double calculateDistance(Parameter X, Parameter Y) {
        double sum = 0.0;
        for (int i = 0; i < mdsDim; i++) {
            double difference = X.getParameterValue(i) - Y.getParameterValue(i);
            sum += difference * difference;
        }
        return Math.sqrt(sum);
    }

    public void makeDirty() {
    }

    // **************************************************************
    // XMLObjectParser
    // **************************************************************

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return ANTIGENIC_TRAIT_LIKELIHOOD;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);

            String fileName = xo.getStringAttribute(FILE_NAME);
            DataTable<double[]> assayTable;
            try {
                assayTable = DataTable.Double.parse(new FileReader(fileName));
            } catch (IOException e) {
                throw new XMLParseException("Unable to read assay data from file, " + fileName);
            }

            // This parameter needs to be linked to the one in the IntegratedMultivariateTreeLikelihood (I suggest that the parameter is created
            // here and then a reference passed to IMTL - which optionally takes the parameter of tip trait values, in which case it listens and
            // updates accordingly.
            CompoundParameter tipTraitParameter = null;
            CompoundParameter virusLocationsParameter = null;
            CompoundParameter serumLocationsParameter = null;

            Parameter mdsPrecision = (Parameter) xo.getElementFirstChild(MDS_PRECISION);

            return new AntigenicTraitLikelihood(treeModel, mdsPrecision, tipTraitParameter, virusLocationsParameter, serumLocationsParameter, assayTable);
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
                new StringAttributeRule(FILE_NAME, "The name of the file containing the assay table"),
                new ElementRule(TreeTraitParserUtilities.TRAIT_PARAMETER, new XMLSyntaxRule[]{
                        new ElementRule(Parameter.class)
                }),
                new ElementRule(MDS_PRECISION, Parameter.class),
                new ElementRule(TreeModel.class),
                new ElementRule(Parameter.class, true)
        };


        public Class getReturnType() {
            return AntigenicTraitLikelihood.class;
        }
    };

    private final double[][] assayTable;

    private final int[] tipIndices;
    private final int[] virusIndices;

    private final CompoundParameter tipTraitParameter;
    private final CompoundParameter virusLocationsParameter;
    private final CompoundParameter serumLocationsParameter;
    private final Parameter mdsParameter;

    private final int totalMeasurementCount;

    // a set of vectors for each virus giving serum indices for which assay data is available
    private final int[][] measuredSerumIndices;

    // a cache of virus to serum distances (serum indices given by array above).
    private final double[][] cachedDistances;

    private final boolean isLeftTruncated;
    private final int mdsDim = 2;
}

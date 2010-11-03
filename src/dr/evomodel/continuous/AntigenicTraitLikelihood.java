package dr.evomodel.continuous;

import dr.evolution.util.Taxon;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.treelikelihood.TreeTraitParserUtilities;
import dr.inference.model.*;
import dr.xml.*;

import java.util.HashMap;
import java.util.Map;

/**
 * @author Andrew Rambaut
 * @author Marc Suchard
 * @version $Id$
 */
public class AntigenicTraitLikelihood extends AbstractModelLikelihood {

    public final static String ANTIGENIC_TRAIT_LIKELIHOOD = "antigenicTraitLikelihood";

    public final static String TIP_TRAIT = "tipTrait";
    public final static String VIRUS_LOCATIONS = "virusLocations";
    public final static String SERUM_LOCATIONS = "serumLocations";

    public AntigenicTraitLikelihood(
            TreeModel tree,
            MultivariateDiffusionModel diffusionModel,
            CompoundParameter tipTraitParameter,
            CompoundParameter virusLocationsParameter,
            CompoundParameter serumLocationsParameter,
            double[][] assayTable,
            String[] virusNames,
            String[] serumNames) {

        super(ANTIGENIC_TRAIT_LIKELIHOOD);

        this.virusNames = virusNames;
        this.serumNames = serumNames;

        int tipCount = tipTraitParameter.getDimension();

        // the total number of viruses is the number of rows in the table
        int virusCount = assayTable.length;
        // the number of sera is the number of columns
        int serumCount = assayTable[0].length;

        if (assayTable.length != virusNames.length) {
            throw new IllegalArgumentException("The number of rows in the assay table doesn't match the virus name list");
        }

        if (assayTable[0].length != serumNames.length) {
            throw new IllegalArgumentException("The number of columns in the assay table doesn't match the serum name list");
        }

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

        for (int i = 0; i < virusCount; i++) {
            virusIndices[i] = -1;

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
                if (!Double.isNaN(assayTable[i][j]) && assayTable[i][j] > 0) {
                    measuredCount ++;
                }
            }

            assayTable[i] = new double[measuredCount];
            measuredSerumIndices[i] = new int[measuredCount];
            cachedDistances[i] = new double[measuredCount];

            int k = 0;
            for (int j = 0; j < serumCount; i++) {
                if (!Double.isNaN(assayTable[i][j]) && assayTable[i][j] > 0) {
                    this.assayTable[i][k] = transform(assayTable[i][k]);
                    measuredSerumIndices[i][k] = j;
                    k ++;
                }
            }
        }

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
        calculateDistances();



        return 0;
    }

    private void calculateDistances() {
        for (int i = 0; i < assayTable.length; i++) {
            for (int j = 0; j < assayTable[i].length; i++) {
                cachedDistances[i][j] = calculateDistance(virusLocationsParameter.getParameter(i), serumLocationsParameter.getParameter(measuredSerumIndices[i][j]));
            }
        }
    }

    private double calculateDistance(Parameter X, Parameter Y) {
        double d1 = X.getParameterValue(0) - Y.getParameterValue(0);
        double d2 = X.getParameterValue(1) - Y.getParameterValue(1);

        return Math.sqrt(d1*d1 + d2*d2);
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

            MultivariateDiffusionModel diffusionModel = (MultivariateDiffusionModel) xo.getChild(MultivariateDiffusionModel.class);
            TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);

            double[][] assayTable = null;
            String[] virusNames = null;
            String[] serumNames = null;

            CompoundParameter tipTraitParameter = null;
            CompoundParameter virusLocationsParameter = null;
            CompoundParameter serumLocationsParameter = null;


            return new AntigenicTraitLikelihood(treeModel, diffusionModel, tipTraitParameter, virusLocationsParameter, serumLocationsParameter, assayTable, virusNames, serumNames);
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
                new StringAttributeRule(TreeTraitParserUtilities.TRAIT_NAME, "The name of the trait for which a likelihood should be calculated"),
                new ElementRule(TreeTraitParserUtilities.TRAIT_PARAMETER, new XMLSyntaxRule[]{
                        new ElementRule(Parameter.class)
                }),
                new ElementRule(MultivariateDiffusionModel.class),
                new ElementRule(TreeModel.class),
                new ElementRule(Parameter.class, true)
        };


        public Class getReturnType() {
            return AntigenicTraitLikelihood.class;
        }
    };

    private final double[][] assayTable;
    private final String[] virusNames;
    private final String[] serumNames;

    private final int[] tipIndices;
    private final int[] virusIndices;

    private final CompoundParameter tipTraitParameter;
    private final CompoundParameter virusLocationsParameter;
    private final CompoundParameter serumLocationsParameter;

    // a set of vectors for each virus giving serum indices for which assay data is available
    private final int[][] measuredSerumIndices;

    // a cache of virus to serum distances (serum indices given by array above).
    private final double[][] cachedDistances;

}

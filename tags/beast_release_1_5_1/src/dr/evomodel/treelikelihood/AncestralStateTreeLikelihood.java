package dr.evomodel.treelikelihood;

import dr.evolution.alignment.PatternList;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.GeneralDataType;
import dr.evolution.tree.NodeAttributeProvider;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.sitemodel.SiteModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Likelihood;
import dr.inference.model.Model;
import dr.math.MathUtils;
import dr.xml.*;

/**
 * @author Marc A. Suchard
 */
public class AncestralStateTreeLikelihood extends TreeLikelihood implements NodeAttributeProvider {

    public static final String RECONSTRUCTING_TREE_LIKELIHOOD = "ancestralTreeLikelihood";
    public static final String RECONSTRUCTION_TAG = "state";
    public static final String TAG_NAME = "tagName";

    private DataType dataType;
    private int[][] reconstructedStates;
    private String tag;


    /**
     * Constructor.
     * Now also takes a DataType so that ancestral states are printed using data codes
     *
     * @param patternList     -
     * @param treeModel       -
     * @param siteModel       -
     * @param branchRateModel -
     * @param useAmbiguities  -
     * @param storePartials   -
     * @param dataType        - need to provide the data-type, so that corrent data characters can be returned
     * @param tag             - string label for reconstruction characters in tree log
     */
    public AncestralStateTreeLikelihood(PatternList patternList, TreeModel treeModel,
                                        SiteModel siteModel, BranchRateModel branchRateModel,
                                        boolean useAmbiguities, boolean storePartials,
                                        DataType dataType,
                                        String tag,
                                        boolean forceRescaling) {
        super(patternList, treeModel, siteModel, branchRateModel, null, useAmbiguities, false, storePartials, false, forceRescaling);
        this.dataType = dataType;
        this.tag = tag;

    }


    public String[] getNodeAttributeLabel() {
        return new String[]{tag};
    }

    public String[] getAttributeForNode(Tree tree, NodeRef node) {
        if (tree != treeModel) {
            throw new RuntimeException("Can only reconstruct states on treeModel given to constructor");
        }

        if (!areStatesRedrawn) {
            redrawAncestralStates();
        }

        return new String[]{formattedState(reconstructedStates[node.getNumber()], dataType)};


    }

    private boolean areStatesRedrawn = false;

    public void redrawAncestralStates() {
        traverseSample(treeModel, treeModel.getRoot(), null);
        areStatesRedrawn = true;
    }

    private boolean checkConditioning = true;


    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == siteModel)
            checkConditioning = true;
        super.handleModelChangedEvent(model, object, index);

    }

    protected double calculateLogLikelihood() {

        if (checkConditioning) {
            final int len = stateCount * stateCount;
            double[] test = new double[len];
            try {
                siteModel.getSubstitutionModel().getTransitionProbabilities(1.0, test);
            } catch (ArithmeticException exception) { // AbstractSubstitutionModel throws numerical errors
                return Double.NEGATIVE_INFINITY;
            }
            for (double d : test) {
                if (d > 1.05) { // ill-conditioned
                    return Double.NEGATIVE_INFINITY;
                }

            }
            checkConditioning = false;
        }

        areStatesRedrawn = false;

        return super.calculateLogLikelihood();
    }

    private static String formattedState(int[] state, DataType dataType) {
        StringBuffer sb = new StringBuffer();
        sb.append("\"");
        if (dataType instanceof GeneralDataType) {
            boolean first = true;
            for (int i : state) {
                if (!first) {
                    sb.append(" ");
                } else {
                    first = false;
                }

                sb.append(dataType.getCode(i));
            }

        } else {
            for (int i : state) {
                sb.append(dataType.getChar(i));
            }
        }
        sb.append("\"");
        return sb.toString();
    }

    /**
     * Traverse (pre-order) the tree sampling the internal node states.
     *
     * @param tree        - TreeModel on which to perform sampling
     * @param node        - current node
     * @param parentState - character state of the parent node to 'node'
     */
    public void traverseSample(TreeModel tree, NodeRef node, int[] parentState) {

        if (reconstructedStates == null)
            reconstructedStates = new int[tree.getNodeCount()][patternCount];

        int nodeNum = node.getNumber();

        NodeRef parent = tree.getParent(node);

        // This function assumes that all partial likelihoods have already been calculated
        // If the node is internal, then sample its state given the state of its parent (pre-order traversal).

        double[] conditionalProbabilities = new double[stateCount];
        int[] state = new int[patternCount];

        if (!tree.isExternal(node)) {

            if (parent == null) {

                double[] rootPartials = getRootPartials();

                // This is the root node
                for (int j = 0; j < patternCount; j++) {

                    System.arraycopy(rootPartials, j * stateCount, conditionalProbabilities, 0, stateCount);
                    try { // TODO This is a hack, need to fix properly
                        state[j] = MathUtils.randomChoicePDF(conditionalProbabilities);
                    } catch (Error e) {
                        System.err.println(e.toString());
                        System.err.println("Please report error to Marc");
                        state[j] = 0;
                    }
                    reconstructedStates[nodeNum][j] = state[j];
                }

            } else {

                // This is an internal node, but not the root
                double[] partialLikelihood = new double[stateCount * patternCount];
                likelihoodCore.getPartials(nodeNum, partialLikelihood);

//				final double branchRate = branchRateModel.getBranchRate(tree, node);
//
//				            // Get the operational time of the branch
//				final double branchTime = branchRate * ( tree.getNodeHeight(parent) - tree.getNodeHeight(node) );
//
//				for (int i = 0; i < categoryCount; i++) {
//
//				                siteModel.getTransitionProbabilitiesForCategory(i, branchTime, probabilities);
//
//				}
//


                if (categoryCount > 1)
                    throw new RuntimeException("Reconstruction not implemented for multiple categories yet.");

                ((AbstractLikelihoodCore) likelihoodCore).getNodeMatrix(nodeNum, 0, probabilities);


                for (int j = 0; j < patternCount; j++) {

                    int parentIndex = parentState[j] * stateCount;
                    int childIndex = j * stateCount;

                    for (int i = 0; i < stateCount; i++)
                        // fixed bug here, index was i, now childIndex + i
                        // is this correct?
                        conditionalProbabilities[i] = partialLikelihood[childIndex + i] * probabilities[parentIndex + i];

                    state[j] = MathUtils.randomChoicePDF(conditionalProbabilities);
                    reconstructedStates[nodeNum][j] = state[j];

                }
            }

            // Traverse down the two child nodes
            NodeRef child1 = tree.getChild(node, 0);
            traverseSample(tree, child1, state);

            NodeRef child2 = tree.getChild(node, 1);
            traverseSample(tree, child2, state);
        } else {

            // This is an external leaf

            ((AbstractLikelihoodCore) likelihoodCore).getNodeStates(nodeNum, reconstructedStates[nodeNum]);

            // Check for ambiguity codes and sample them

            for (int j = 0; j < patternCount; j++) {

                final int thisState = reconstructedStates[nodeNum][j];

                if (dataType.isAmbiguousState(thisState)) {

                    int parentIndex = parentState[j] * stateCount;
                    ((AbstractLikelihoodCore) likelihoodCore).getNodeMatrix(nodeNum, 0, probabilities);
                    System.arraycopy(probabilities, parentIndex, conditionalProbabilities, 0, stateCount);
                    reconstructedStates[nodeNum][j] = MathUtils.randomChoicePDF(conditionalProbabilities);

                }

            }
        }
    }


    /**
     * The XML parser
     */
    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String[] getParserNames() {
            return new String[] {
                    getParserName(),"beast_"+getParserName()
            };
        }

        public String getParserName() {
            return RECONSTRUCTING_TREE_LIKELIHOOD;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            boolean useAmbiguities = xo.getAttribute(USE_AMBIGUITIES, false);
            boolean storePartials = xo.getAttribute(STORE_PARTIALS, true);

            PatternList patternList = (PatternList) xo.getChild(PatternList.class);
            TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
            SiteModel siteModel = (SiteModel) xo.getChild(SiteModel.class);

            BranchRateModel branchRateModel = (BranchRateModel) xo.getChild(BranchRateModel.class);

            DataType dataType = ((SubstitutionModel) xo.getChild(SubstitutionModel.class)).getDataType();

            // default tag is RECONSTRUCTION_TAG
            String tag = xo.getAttribute(TAG_NAME, RECONSTRUCTION_TAG);

            boolean forceRescaling = xo.getAttribute(TreeLikelihood.FORCE_RESCALING,false);

            return new AncestralStateTreeLikelihood(patternList, treeModel, siteModel,
                    branchRateModel, useAmbiguities, storePartials, dataType, tag, forceRescaling);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element represents the likelihood of a patternlist on a tree given the site model.";
        }

        public Class getReturnType() {
            return Likelihood.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                AttributeRule.newBooleanRule(USE_AMBIGUITIES, true),
                AttributeRule.newBooleanRule(STORE_PARTIALS, true),
                AttributeRule.newStringRule(TAG_NAME, true),
                AttributeRule.newBooleanRule(TreeLikelihood.FORCE_RESCALING,true),
                new ElementRule(PatternList.class),
                new ElementRule(TreeModel.class),
                new ElementRule(SiteModel.class),
                new ElementRule(BranchRateModel.class, true),
                new ElementRule(SubstitutionModel.class)
        };
    };

}

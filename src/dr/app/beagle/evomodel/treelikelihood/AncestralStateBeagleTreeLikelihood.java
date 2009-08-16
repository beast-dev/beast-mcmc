package dr.app.beagle.evomodel.treelikelihood;

import dr.evolution.alignment.PatternList;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.GeneralDataType;
import dr.evolution.tree.NodeAttributeProvider;
import dr.evolution.tree.Tree;
import dr.evolution.tree.NodeRef;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.treelikelihood.AncestralStateTreeLikelihood;
import dr.evomodel.sitemodel.SiteModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.app.beagle.evomodel.sitemodel.BranchSiteModel;
import dr.app.beagle.evomodel.sitemodel.SiteRateModel;
import dr.math.MathUtils;
import beagle.Beagle;

/**
 * @author Marc Suchard
 * @author Andrew Rambaut
 */

public class AncestralStateBeagleTreeLikelihood extends BeagleTreeLikelihood implements NodeAttributeProvider {

    public AncestralStateBeagleTreeLikelihood(PatternList patternList, TreeModel treeModel,
                                              BranchSiteModel branchSiteModel, SiteRateModel siteRateModel,
                                              BranchRateModel branchRateModel, boolean useAmbiguities,
                                              boolean alwaysRescale,
                                               DataType dataType,
                                            String tag) {

        super(patternList, treeModel, branchSiteModel, siteRateModel, branchRateModel, useAmbiguities, alwaysRescale);

        this.dataType = dataType;
        this.tag = tag;

//        branchSiteModel.


//        oldSiteModel = new SiteModel();

//        oldAncestralLikelihood = new AncestralStateTreeLikelihood()

//                                                        SiteModel siteModel, BranchRateModel branchRateModel,
//                                        boolean useAmbiguities, boolean storePartials,
//                                        DataType dataType,
//                                        String tag) {

        probabilities = new double[stateCount*stateCount];
        partials = new double[stateCount*patternCount];
        rootPartials = new double[stateCount*patternCount];
        cumulativeScaleBuffers = new int[nodeCount][];
        scaleBufferIndex = getScaleBufferCount() - 1;
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

    @Override
     protected int getScaleBufferCount() {
        return internalNodeCount + 2;
    }

     private boolean areStatesRedrawn = false;

     public void redrawAncestralStates() {
         // Setup cumulate scale buffers
         traverseCollectScaleBuffers(treeModel, treeModel.getRoot());
         // Sample states
         traverseSample(treeModel, treeModel.getRoot(), null);
         areStatesRedrawn = true;
     }

    protected double calculateLogLikelihood() {
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

//    private double[] getRootPartials() {
//        int rootNum = treeModel.getRoot().getNumber();
//        beagle.resetScaleFactors(scaleBufferIndex);
//        beagle.accumulateScaleFactors(cumulativeScaleBuffers[rootNum],cumulativeScaleBuffers[rootNum].length,scaleBufferIndex);
//        beagle.getPartials(partialBufferHelper.getOffsetIndex(rootNum),scaleBufferIndex,rootPartials);
//        return rootPartials;
//    }

    private void getPartials(int number, double[] partials) {
        int cumulativeBufferIndex = Beagle.NONE;
        if (useScaleFactors) {
            cumulativeBufferIndex = scaleBufferIndex;
            beagle.resetScaleFactors(cumulativeBufferIndex);
            beagle.accumulateScaleFactors(cumulativeScaleBuffers[number],cumulativeScaleBuffers[number].length,cumulativeBufferIndex);
        }
        beagle.getPartials(partialBufferHelper.getOffsetIndex(number),cumulativeBufferIndex,partials);
    }

    private void getMatrix(int matrixNum, double[] probabilities) {
//        beagle.getMatrix(); // TODO Need to add accessor to BEAGLE API.  Oh no!!!
    }

    private void getStates(int tipNum, int[] states)  {
        // TODO
    }

    public int traverseCollectScaleBuffers(TreeModel tree, NodeRef node) {

        if (true) // Currently do nothing
            return 0;

        if (tree.isExternal(node))
            return 0;

        int nodeNum = node.getNumber();

        NodeRef child0 = tree.getChild(node,0);
        NodeRef child1 = tree.getChild(node,1);

        int len0 = traverseCollectScaleBuffers(tree,child0);
        int len1 = traverseCollectScaleBuffers(tree,child1);
        int thisLen = len0 + len1 + 1;
        int offset = 0;

        int[] scaleBuffer = new int[thisLen];
        if (len0 > 0) {
            System.arraycopy(cumulativeScaleBuffers[child0.getNumber()],0,scaleBuffer,offset,len0);
            offset += len0;
        }
        if (len1 > 0) {
            System.arraycopy(cumulativeScaleBuffers[child1.getNumber()],0,scaleBuffer,offset,len1);
            offset += len1;
        }
        scaleBuffer[offset] = scaleBufferHelper.getOffsetIndex(nodeNum - tipCount);        
        cumulativeScaleBuffers[nodeNum] = scaleBuffer;

        return thisLen;
    }


    public void traverseSample(TreeModel tree, NodeRef node, int[] parentState) {

        if (reconstructedStates == null)
            reconstructedStates = new int[tree.getNodeCount()][patternCount];

        if (true)  // Currently do nothing
            return;
        
        int nodeNum = node.getNumber();

        NodeRef parent = tree.getParent(node);

        // This function assumes that all partial likelihoods have already been calculated
        // If the node is internal, then sample its state given the state of its parent (pre-order traversal).

        double[] conditionalProbabilities = new double[stateCount];
        int[] state = new int[patternCount];

        if (!tree.isExternal(node)) {

            if (parent == null) {

                 // This is the root node
                getPartials(nodeNum,partials);
                for (int j = 0; j < patternCount; j++) {

                    System.arraycopy(partials, j * stateCount, conditionalProbabilities, 0, stateCount);
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
                getPartials(nodeNum,partialLikelihood);  // TODO fix number

                if (categoryCount > 1)
                    throw new RuntimeException("Reconstruction not implemented for multiple categories yet.");

                getMatrix(nodeNum*categoryCount,probabilities); // TODO fix number

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

            getStates(nodeNum, reconstructedStates[nodeNum]); // TODO fix number

            // Check for ambiguity codes and sample them

            for (int j = 0; j < patternCount; j++) {

                final int thisState = reconstructedStates[nodeNum][j];

                if (dataType.isAmbiguousState(thisState)) {

                    int parentIndex = parentState[j] * stateCount;
                    getMatrix(nodeNum*categoryCount,probabilities); // TODO fix number
                    System.arraycopy(probabilities, parentIndex, conditionalProbabilities, 0, stateCount);
                    reconstructedStates[nodeNum][j] = MathUtils.randomChoicePDF(conditionalProbabilities);
                }
            }
        }
    }

    private DataType dataType;
    private int[][] reconstructedStates;
    private String tag;

    private AncestralStateTreeLikelihood oldAncestralLikelihood;
    private SiteModel oldSiteModel;
    private SubstitutionModel oldSubstitutionModel;

    private double[] probabilities;
    private double[] partials;
    private double[] rootPartials;
    private int[][] cumulativeScaleBuffers;
    private int scaleBufferIndex;
}

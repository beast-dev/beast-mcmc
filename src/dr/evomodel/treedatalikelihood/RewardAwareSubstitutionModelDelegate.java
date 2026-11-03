package dr.evomodel.treedatalikelihood;

import beagle.Beagle;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.branchmodel.TransitionMatrixProviderBranchModel;

/**
 * RewardAwareSubstitutionModelDelegate
 *
 * Contract:
 *  - We do NOT use eigen decompositions or BEAGLE's eigen-based updateTransitionMatrices.
 *  - Instead, the BranchModel provides a full transition matrix W per branch (NodeRef),
 *    and we push W directly into BEAGLE transition-matrix buffers.
 *
 * Critical indexing rules:
 *  1) BEAGLE transition matrices are stored in buffers indexed by the delegate's matrix buffer helper:
 *         beagleMatrixIndex = getMatrixIndex(branchIndex)
 *  2) When flipBuffers == true, we must flip the transition-matrix buffers for those branch indices
 *     before writing, so we write into the "current" side that BEAST will read.
 */

/*
 @author Filippo Monti
 */

public final class RewardAwareSubstitutionModelDelegate extends SubstitutionModelDelegate {

    private final TransitionMatrixProviderBranchModel branchModel;
    private final Tree tree;

    public RewardAwareSubstitutionModelDelegate(Tree tree,
                                                TransitionMatrixProviderBranchModel branchModel,
                                                int partitionNumber,
                                                int bufferPoolSize,
                                                PreOrderSettings settings) {
        super(tree, branchModel, partitionNumber, bufferPoolSize, settings);
        this.branchModel = branchModel;
        this.tree = tree;
    }

    @Override
    public void updateSubstitutionModels(Beagle beagle, boolean flipBuffers) {
        // Intentionally no-op: transition matrices are provided directly by branchModel.
    }

    /**
     * Push transition matrices directly into BEAGLE buffers.
     *
     * @param branchIndices node-indices (compatible with operations arrays) whose matrices need updating
     * @param edgeLength ignored (branchModel encodes all time/reward dependence internally)
     * @param updateCount number of items in branchIndices/edgeLength to process
     * @param flipBuffers whether BEAST is requesting a flip of the matrix buffer offsets
     */
    @Override
    public void updateTransitionMatrices(Beagle beagle,
                                         int[] branchIndices,
                                         double[] edgeLength,
                                         int updateCount,
                                         boolean flipBuffers) {

        // Respect the parent's double-buffering scheme.
        if (flipBuffers) {
            flipTransitionMatrices(branchIndices, updateCount);
        }

        for (int i = 0; i < updateCount; i++) {
            final int branchIndex = branchIndices[i];

            // branchIndex is a node-index into tree.getNode(branchIndex)
            final NodeRef node = tree.getNode(branchIndex);

            // Full transition matrix for this branch (row-major, length = dim*dim)
            final double[] W = branchModel.getTransitionMatrix(node);

            // write into the delegate-managed matrix buffer index
            final int beagleMatrixIndex = getMatrixIndex(branchIndex);
            beagle.setTransitionMatrix(beagleMatrixIndex, W, 0);
        }
    }

    /**
     * debugging utilities: for sanity checking during debugging.
     * TODO check node index != node.getNumber().
     */
    @SuppressWarnings("unused")
    public static void printNodeIndexVsNumber(Tree tree) {
        for (int i = 0; i < tree.getNodeCount(); i++) {
            NodeRef node = tree.getNode(i);
            System.out.println("treeIndex=" + i + " node.getNumber()=" + node.getNumber()
                    + (tree.isRoot(node) ? " (root)" : ""));
        }
    }

    public static void printBranchOrder(Tree tree) {
        for (int i = 0; i < tree.getNodeCount(); i++) {
            NodeRef node = tree.getNode(i);
            if (tree.isRoot(node)) continue;

            int childNr = node.getNumber();
            NodeRef parent = tree.getParent(node);
            int parentNr = parent.getNumber();

            String childName = tree.isExternal(node) ? tree.getNodeTaxon(node).getId() : ("internal#" + childNr);
            String parentName = tree.isExternal(parent) ? tree.getNodeTaxon(parent).getId() : ("internal#" + parentNr);

            double t = tree.getBranchLength(node);

            System.out.println(
                    "nodeNr=" + childNr +
                            " parentNr=" + parentNr +
                            " child=" + childName +
                            " parent=" + parentName +
                            " t=" + t
            );
        }
    }


}

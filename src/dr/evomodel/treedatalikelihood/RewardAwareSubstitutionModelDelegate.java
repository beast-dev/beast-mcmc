package dr.evomodel.treedatalikelihood;

import beagle.Beagle;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.branchmodel.TransitionMatrixProviderBranchModel;
import dr.evomodel.substmodel.SubstitutionModel;

public class RewardAwareSubstitutionModelDelegate extends SubstitutionModelDelegate {
    private final TransitionMatrixProviderBranchModel branchModel;
    private final Tree tree;

    public RewardAwareSubstitutionModelDelegate(Tree tree, TransitionMatrixProviderBranchModel branchModel,
                                                int partitionNumber, int bufferPoolSize,
                                                PreOrderSettings settings) {
        super(tree, branchModel, partitionNumber, bufferPoolSize, settings);
        this.branchModel = branchModel;
        this.tree = tree;
    }

    @Override
    public void updateSubstitutionModels(Beagle beagle, boolean flipBuffers) {
        //   do nothing
    }

    @Override
    public void updateTransitionMatrices(Beagle beagle, int[] branchIndices, double[] edgeLength, int updateCount, boolean flipBuffers) {

        for (int i = 0; i < updateCount; i++) {
            NodeRef node = tree.getNode(branchIndices[i]);
            final double[] W = branchModel.getTransitionMatrix(node);
            beagle.setTransitionMatrix(branchIndices[i], W, 0); // TODO check this indexing
        }

    }// END: updateTransitionMatrices

    // TODO this was same
    @Override
    public int getSubstitutionModelCount() {
        return tree.getNodeCount() - 1;
    }

    //TODO THIS IS SIMILAR
    @Override
    public SubstitutionModel getSubstitutionModel(int index) {
        return null;
    }


}

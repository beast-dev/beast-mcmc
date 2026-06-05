package dr.evomodel.branchmodel;

import dr.evolution.tree.NodeRef;

public interface TransitionMatrixProviderBranchModel extends BranchModel {

        /**
        * Returns the transition matrix for the given branch.
        *
        * @param branch the branch
        * @return the transition matrix
        */
        double[] getTransitionMatrix(final NodeRef branch);
}

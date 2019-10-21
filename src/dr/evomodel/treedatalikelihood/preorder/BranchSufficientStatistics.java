package dr.evomodel.treedatalikelihood.preorder;

import dr.math.matrixAlgebra.missingData.PermutationIndices;

/**
 * @author Marc A. Suchard
 */
public class BranchSufficientStatistics {

    private final NormalSufficientStatistics below;
    private final MatrixSufficientStatistics branch;
    private final NormalSufficientStatistics above;

    BranchSufficientStatistics(NormalSufficientStatistics below,
                               MatrixSufficientStatistics branch,
                               NormalSufficientStatistics above) {
        this.below = below;
        this.branch = branch;
        this.above = above;
    }

    public NormalSufficientStatistics getBelow() { return below; }

    public MatrixSufficientStatistics getBranch() { return branch; }

    public NormalSufficientStatistics getAbove() { return above; }

    public String toString() { return below + " / " + branch + " / " + above; }

    public String toVectorizedString() {
        return below.toVectorizedString() + " / " + branch.toVectorizedString() + " / " + above.toVectorizedString();
    }

    public int[] getMissing() {
        PermutationIndices indices = new PermutationIndices(getBelow().getRawPrecision());
        return indices.getZeroIndices();
    }
}

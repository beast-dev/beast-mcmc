package dr.evomodel.treedatalikelihood.preorder;

/**
 * @author Marc A. Suchard
 */
public class BranchSufficientStatistics {

    private final NormalSufficientStatistics child;
    private final MatrixSufficientStatistics branch;
    private final NormalSufficientStatistics parent;

    BranchSufficientStatistics(NormalSufficientStatistics child,
                               MatrixSufficientStatistics branch,
                               NormalSufficientStatistics parent) {
        this.child = child;
        this.branch = branch;
        this.parent = parent;
    }

    public NormalSufficientStatistics getChild() { return child; }

    public MatrixSufficientStatistics getBranch() { return branch; }

    public NormalSufficientStatistics getParent() { return parent; }

    public String toString() { return child + " / " + branch + " / " + parent; }

    public String toVectorizedString() {
        return child.toVectorizedString() + " / " + branch.toVectorizedString() + " / " + parent.toVectorizedString();
    }
}

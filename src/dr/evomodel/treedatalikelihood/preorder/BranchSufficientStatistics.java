package dr.evomodel.treedatalikelihood.preorder;

/**
 * @author Marc A. Suchard
 */
public class BranchSufficientStatistics {

    private final NormalSufficientStatistics child;
    private final NormalSufficientStatistics branch;
    private final NormalSufficientStatistics parent;

    BranchSufficientStatistics(NormalSufficientStatistics child,
                               NormalSufficientStatistics branch,
                               NormalSufficientStatistics parent) {
        this.child = child;
        this.branch = branch;
        this.parent = parent;
    }

    public NormalSufficientStatistics getChild() { return child; }

    public NormalSufficientStatistics getBranch() { return branch; }

    public NormalSufficientStatistics getParent() { return parent; }

    public String toString() { return child + " / " + branch + " / " + parent; }

    public String toVectorizedString() {
        return child.toVectorizedString() + " / " + branch.toVectorizedString() + " / " + parent.toVectorizedString();
    }
}

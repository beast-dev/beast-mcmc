package dr.evomodel.coalescent;
import dr.evomodel.bigfasttree.BigFastTreeIntervals;

// this builds a timeline with coalescent and sampling events for a single tree
public abstract class SingleTreeTimeline {
    protected final BigFastTreeIntervals treeIntervals;
    final int nNodes;

    public SingleTreeTimeline(BigFastTreeIntervals treeIntervals) {
        this.treeIntervals = treeIntervals;
        this.nNodes = treeIntervals.getIntervalCount() + 1;
    }

    protected int getnNodes() { return nNodes; }
    protected BigFastTreeIntervals getTreeIntervals() { return treeIntervals; }
}
package dr.evomodel.coalescent;

// this builds a timeline with coalescent and sampling events for a single tree
public abstract class SingleTreeTimeline {
    protected final TreeIntervals treeIntervals;
    final int nNodes;

    public SingleTreeTimeline(TreeIntervals treeIntervals) {
        this.treeIntervals = treeIntervals;
        this.nNodes = treeIntervals.getIntervalCount() + 1;
    }

    protected int getnNodes() { return nNodes; }
    protected TreeIntervals getTreeIntervals() { return treeIntervals; }
}
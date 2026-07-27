package dr.evomodelxml.coalescent;
import dr.evomodel.bigfasttree.BigFastTreeIntervals;

// this builds a timeline with coalescent and sampling events for a single tree
public class SingleTreeNodesTimeline {
    private final int nNodes;
    private final BigFastTreeIntervals treeIntervals;

    private final double[] timeLine;
    private final boolean[] flagCoalescentEvent;
    private int[] numLineages;

    public SingleTreeNodesTimeline(BigFastTreeIntervals treeIntervals) {
        this.treeIntervals = treeIntervals;
        this.nNodes = treeIntervals.getIntervalCount() + 1;
        this.timeLine = new double[nNodes];
        this.flagCoalescentEvent = new boolean[nNodes];
        makeLine();
    }

    private void makeLine() {
        timeLine[0] = treeIntervals.getStartTime();
        flagCoalescentEvent[0] = false;
        for (int nodeIndex = 1; nodeIndex < nNodes; nodeIndex++) {
            timeLine[nodeIndex] = treeIntervals.getIntervalTime(nodeIndex);
            flagCoalescentEvent[nodeIndex] =
                    String.valueOf(treeIntervals.getIntervalType(nodeIndex - 1)).equals("coalescent"); //TODO this is hard coded ...
        }
    }

    protected int getnNodes() { return nNodes; }
    protected boolean[] getFlagCoalescentEvent() { return flagCoalescentEvent; }
    protected double[] getTimeLine() { return timeLine; }
}
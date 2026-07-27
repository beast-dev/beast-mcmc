package dr.evomodel.coalescent;
import dr.evomodel.bigfasttree.BigFastTreeIntervals;

// this builds a timeline with coalescent and sampling events for a single tree
public class SingleTreeNodesTimeline extends SingleTreeTimeline {
    private final int nNodes;

    private final double[] timeLine;
    private final boolean[] flagCoalescentEvent;
    private int[] numLineages;

    public SingleTreeNodesTimeline(BigFastTreeIntervals treeIntervals) {
        super(treeIntervals);
        this.nNodes = getnNodes();
        this.timeLine = new double[nNodes];
        this.numLineages = new int[nNodes];
        this.flagCoalescentEvent = new boolean[nNodes];
        makeLine();
    }

    protected void makeLine() {
        timeLine[0] = treeIntervals.getStartTime();
        flagCoalescentEvent[0] = false;
        for (int nodeIndex = 1; nodeIndex < nNodes; nodeIndex++) {
            timeLine[nodeIndex] = treeIntervals.getIntervalTime(nodeIndex);
            numLineages[nodeIndex] = treeIntervals.getLineageCount(nodeIndex);
            flagCoalescentEvent[nodeIndex] =
                    String.valueOf(treeIntervals.getIntervalType(nodeIndex - 1)).equals("coalescent"); //TODO this is hard coded ...
        }
    }


    protected int getnNodes() { return nNodes; }
    protected boolean[] getFlagCoalescentEvent() { return flagCoalescentEvent; }
    protected double[] getTimeLine() { return timeLine; }
}
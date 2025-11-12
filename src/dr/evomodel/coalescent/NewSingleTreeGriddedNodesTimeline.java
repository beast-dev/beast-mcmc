package dr.evomodel.coalescent;

import dr.inference.model.Parameter;

import java.util.Arrays;

// this builds a timeline with coalescent and sampling events for a single tree with grid points
public class NewSingleTreeGriddedNodesTimeline extends NewSingleTreeTimeline {
    private final int nNodes;
    private final double[] timeLine;
    private final boolean[] flagCoalescentEvent;

    private final Parameter gridPoints;
    private int[] gridIndices;

    private double[] mergedTimeLine;
    private int[] mergedNumLineages;
    private int[] numCoalEvents;

    public NewSingleTreeGriddedNodesTimeline(TreeIntervals treeIntervals, Parameter gridPoints) {
        super(treeIntervals);
        this.nNodes = getnNodes();
        this.gridPoints = gridPoints;
        this.timeLine = new double[nNodes];
        this.flagCoalescentEvent = new boolean[nNodes];
        makeNodesLine();

        if (gridPoints != null) { // "merged" =  nodes' times and grid points
            this.mergedTimeLine = new double[nNodes + gridPoints.getDimension()];
            this.gridIndices = new int[gridPoints.getDimension()];
            numCoalEvents = new int[gridPoints.getDimension() + 1]; // "+1" to account for the events after the last grid point
            mergedNumLineages = new int[nNodes + gridPoints.getDimension() + 1];
            integrateGridPoints(gridPoints.getParameterValues());
        }

    }

    protected void makeNodesLine() {
        timeLine[0] = treeIntervals.getStartTime();
        flagCoalescentEvent[0] = false;
        for (int nodeIndex = 1; nodeIndex < nNodes; nodeIndex++) {
            timeLine[nodeIndex] = treeIntervals.getIntervalTime(nodeIndex);
            flagCoalescentEvent[nodeIndex] =
                    String.valueOf(treeIntervals.getIntervalType(nodeIndex - 1)).equals("coalescent"); //TODO this is hard coded ...
        }
    }

    private void integrateGridPoints(double[] gridPointsVector) {
        int nGridPoints = gridPointsVector.length;
        Arrays.fill(mergedNumLineages, 0);
        Arrays.fill(numCoalEvents, 0);

        int i = 0, j = 0;
        while (i < nNodes && j < nGridPoints) {
            if (timeLine[i] <= gridPointsVector[j]) { // grid points are set after node times if equal
                mergedTimeLine[i + j] = timeLine[i];
                computeNumLineages(i + j, flagCoalescentEvent[i]);
                if (flagCoalescentEvent[i]) numCoalEvents[j] += 1;
                i++;
            } else {
                mergedTimeLine[i + j] = gridPointsVector[j];
                gridIndices[j] = i + j;
                mergedNumLineages[i + j] = mergedNumLineages[i + j - 1];
                j++;
            }
        }

        // Add remaining elements
        while (i < nNodes) {
            mergedTimeLine[i + j] = timeLine[i];
            computeNumLineages(i + j, flagCoalescentEvent[i]);
            if (flagCoalescentEvent[i]) numCoalEvents[j] += 1;
            i++;
        }
        while (j < nGridPoints) {
            mergedTimeLine[i + j] = gridPointsVector[j];
            gridIndices[j] = i + j;
            mergedNumLineages[i + j] = mergedNumLineages[i + j - 1];
            j++;
        }
    }

    private void computeNumLineages(int index, boolean isCoalescentEvent) {
        if (isCoalescentEvent) {
            mergedNumLineages[index] = mergedNumLineages[index - 1] - 1;
        } else {
            if (index == 0) { // the first node is a tip node
                mergedNumLineages[index] = 1;
            } else {
                mergedNumLineages[index] = mergedNumLineages[index - 1] + 1;
            }
        }
    }

    public double[] getMergedTimeLine() { return mergedTimeLine; }
    public int[] getMergedNumLineages() {return mergedNumLineages;}
    public int[] getNumCoalEvents() {return numCoalEvents;}
    public int[] getGridIndices() {return gridIndices;}
}
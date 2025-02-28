package dr.evomodel.branchratemodel;

import dr.evomodel.bigfasttree.BigFastTreeIntervals;
import dr.evomodel.coalescent.SingleTreeTimeline;
import dr.inference.model.Parameter;

// this builds a timeline with coalescent and sampling events for a single tree
public class NPRatesTimeline extends SingleTreeTimeline {
    private final int nNodes;
    private final Parameter gridPoints;
    private final int nGridPoints;

    private final double[] nodesTimeLine; // only with nodes
    private final double[] fullTimeLine; // nodes and grid points

//    the pair of the following two vector creates a mapping child -> parent
//    private int[] childNodes; // lists all child nodes
    private int[] parentNodes; // lists all parent nodes

    private int[] gridAndChildNodesFlags; // flags grid points
    private double[][] branchesIntersections; // for each interval induced by the grid points store the sub-intervals the intersect each branch

    public NPRatesTimeline(BigFastTreeIntervals treeIntervals, Parameter gridPoints) {
        super(treeIntervals);
        this.nNodes = getnNodes();

        this.gridPoints = gridPoints;
        this.nGridPoints = gridPoints.getDimension();

        this.nodesTimeLine = new double[nNodes];

        this.fullTimeLine = new double[nNodes + nGridPoints];
        this.gridAndChildNodesFlags = new int[nNodes + nGridPoints];

//        this.childNodes = new int[nNodes];
        this.parentNodes = new int[nNodes]; // TODO correct dimension and update code below

        this.branchesIntersections = new double[nGridPoints + 1][nNodes - 1]; // for each interval induced by the grid points store the sub-intervals the intersect each branch

        makeNodesTimeLine(); //  create the nodesTimeLine and the mapping "childNodes to parentNodes"
        makeFullTimeLine(); // second I create the fullTimeLine considering both nodes and grid points
        makeGridNodesIntersections(); // make the mapping between intervals induced by grid points and branches
    }
// TODO this code is under the assumption that the branches are named based on the child node
    private void makeNodesTimeLine() {
        nodesTimeLine[0] = treeIntervals.getStartTime();
        for (int nodeIndex = 1; nodeIndex < nNodes; nodeIndex++) {
            nodesTimeLine[nodeIndex] = treeIntervals.getIntervalTime(nodeIndex);
//            childNodes[nodeIndex] = 1; //TODO get the number associated with the branch
            parentNodes[nodeIndex] = 0; //TODO get the parent of the previous node--> should do from treeModel??
        }
    }

    private void makeFullTimeLine() {
        int i = 0, j = 0;
        while (i < nNodes && j < nGridPoints) {
            if (gridPoints.getParameterValue(j) <= nodesTimeLine[i]) { // grid points put before nodes
                fullTimeLine[i + j] = gridPoints.getParameterValue(j);
                gridAndChildNodesFlags[i + j] = - j; // negative values for grid points
                j++;
            } else {
                fullTimeLine[i + j] = nodesTimeLine[i];
                // TODO this is under the assumpting that nodes are named in order of their height
                gridAndChildNodesFlags[i + j] = i; // positive values for nodes

                i++;
            }
        }
        while (i < nNodes) {
            fullTimeLine[i + j] = nodesTimeLine[i];
            gridAndChildNodesFlags[i + j] = i;
            i++;
        }
        while (j < nGridPoints) {
            fullTimeLine[i + j] = gridPoints.getParameterValue(j);
            gridAndChildNodesFlags[i + j] = - j; // negative values for grid points
            j++;
        }
    }

    private void makeGridNodesIntersections() {
        int index = 0;
        while (index < fullTimeLine.length) {
            if (gridAndChildNodesFlags[index] > 0) { // if it is a (child) node, and so indexes a branch
                double t = fullTimeLine[index];
                int k = index + 1; // jump to the next index
                // now compute the intersection with the next grid point
                // and skip if it is a node that is not the child node's parent
                while (gridAndChildNodesFlags[k] < 0 || gridAndChildNodesFlags[k] != parentNodes[gridAndChildNodesFlags[index]]) { // compute intersection if the index is a grid point
                    if (gridAndChildNodesFlags[k] < 0) {
                        branchesIntersections[-gridAndChildNodesFlags[k]][index] = fullTimeLine[k] - t;
                        t = fullTimeLine[k];
                    }
                    k++;
                }
                double tParent = fullTimeLine[k];
                while (gridAndChildNodesFlags[k] > 0) k++; // get to the next grid point
                if (k == fullTimeLine.length) { // check if going after the last grid point
                    branchesIntersections[nGridPoints + 1][index] = tParent - t;
                } else {
                    branchesIntersections[-gridAndChildNodesFlags[k]][index] = tParent - t;
                }

            }
        }
    }

    protected int getnNodes() { return nNodes; }
    protected double[] getNodesTimeLine() { return nodesTimeLine; }
    protected double[] getFullTimeLine() { return fullTimeLine; }
    protected int[] getGridAndChildNodesFlags() { return gridAndChildNodesFlags; }
    protected double[][] getBranchesIntersections() { return branchesIntersections; }
}

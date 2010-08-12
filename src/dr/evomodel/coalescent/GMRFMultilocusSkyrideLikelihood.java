package dr.evomodel.coalescent;

import dr.evolution.coalescent.TreeIntervals;
import dr.evolution.tree.Tree;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.MatrixParameter;
import dr.inference.model.Model;
import dr.inference.model.Parameter;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Marc A. Suchard
 * @author Philippe Lemey
 * @author Erik W. Bloomquist
 */

public class GMRFMultilocusSkyrideLikelihood extends GMRFSkyrideLikelihood implements MultiLociTreeSet {

    	public GMRFMultilocusSkyrideLikelihood(List<Tree> trees,
                                               Parameter popParameter,
                                               Parameter groupParameter,
                                               Parameter precParameter,
	                                           Parameter lambda,
                                               Parameter beta,
                                               MatrixParameter dMatrix,
	                                           boolean timeAwareSmoothing) {
            super(trees, popParameter, groupParameter, precParameter, lambda, beta, dMatrix, timeAwareSmoothing);
        }

    protected void setTree(List<Tree> treeList) {
        treesSet = this;
        this.treeList = treeList;
        intervalsList = new ArrayList<TreeIntervals>();
        for (Tree tree : treeList) {
            intervalsList.add(new TreeIntervals(tree));
            if (tree instanceof TreeModel) {
                addModel((TreeModel) tree);
            }
        }
    }

    protected int getCorrectFieldLength() {
        // TODO Fix!
        int tips = 0;
        for (Tree tree : treeList) {
            tips += tree.getExternalNodeCount();
        }
        return tips - treeList.size();
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {

        if (model instanceof TreeModel) {
            TreeModel treeModel = (TreeModel) model;
            int tn = treeList.indexOf(treeModel);
            if (tn >= 0) {
                intervalsList.get(tn).setIntervalsUnknown();
            }
        }
        
        super.handleModelChangedEvent(model, object, index);
    }

    public void initializationReport() {
		System.out.println("Creating a GMRF smoothed skyride model for multiple loci:");
		System.out.println("\tPopulation sizes: " + popSizeParameter.getDimension());
		System.out.println("\tIf you publish results using this model, please reference: ");
        System.out.println("\t\tMinin, Bloomquist and Suchard (2008) Molecular Biology and Evolution, 25, 1459-1471, and");
        System.out.println("\t\tSuchard, Drummond and Lemey (in preparation).");
	}

    public void wrapSetupIntervals() {
        // Do nothing
    }

    protected void setupSufficientStatistics() {


//What we want:

//-Matrix with number of lineages present in each tree during each interval
//-length of each interval
//-vector which has -1 if corresponding interval does not end with coalescent event, and otherwise has i, where i is the tree in which the coalescent event took place


//Add the following method to TreeIntervals.java:

//get time of node which marks start of interval
//public double getIntervalTime(int i){
//	return times[indices[i]];
//}

        int numTrees = treeList.size();

        int N = 0;

        for (int i = 0; i < numTrees; i++) {
            N = N + intervalsList.get(i).getIntervalCount();
        }


        int[][] lineageCounts = new int[numTrees][N];
        double[] intervalLengths = new double[N];
        int[] coalEvent = new int[N];

        //vector keeping track of which intervals are next in each tree
        //should be initialized to all zeros
        int nextIntervals[] = new int[numTrees];

        int treeFinished[] = new int[numTrees]; //vector indicating whether we have reached end of tree
        int numTreesFinished = 0;


        int index = 0;

        int currentTree = 0;
        double currentTime = intervalsList.get(0).getIntervalTime(0);
        double propCurrent;

        for (int i = 1; i < numTrees; i++) {
            propCurrent = intervalsList.get(i).getIntervalTime(0);
            if (propCurrent < currentTime) {
                currentTree = i;
                currentTime = propCurrent;
            }
        }


        for (TreeIntervals treeInterval : intervalsList) {
            System.err.println(treeInterval);
        }


        nextIntervals[currentTree] = nextIntervals[currentTree] + 1;

        System.err.println("Start of while-loop");



        while (numTreesFinished < numTrees) {

            int nextTree = 0;
            double nextTime = 0;
            double length = 0;

            while (nextTree < numTrees) {
                if (treeFinished[nextTree] == 0) {
                    nextTime = intervalsList.get(nextTree).getIntervalTime(nextIntervals[nextTree]);
                    if (nextTime > currentTime) {  // to be sure there are no intervals of length zero
                        length = nextTime - currentTime;
                        break;
                    }
                } else {
                    nextTree++;
                }
            }


            double propNext;
            int startIndex = nextTree + 1;   //for efficiency

            for (int i = startIndex; i < numTrees; i++) {
                if (treeFinished[i] == 0) {
                    propNext = intervalsList.get(i).getIntervalTime(nextIntervals[i]);
                    if ((propNext - currentTime) < length) {
                        nextTree = i;
                        nextTime = propNext;
                        length = nextTime - currentTime;
                    }
                }
            }


            int tempNode;

            for (int i = 0; i < numTrees; i++) {
                if (treeFinished[i] == 0) {
                    tempNode = nextIntervals[i] - 1;
                    if (tempNode < 0) {
                        lineageCounts[i][index] = 0;
                    } else {
                        lineageCounts[i][index] = intervalsList.get(i).getLineageCount(tempNode);
                    }
                }
            }


//what if coalescent events in different trees happen at the same time?
            if (intervalsList.get(nextTree).getCoalescentEvents(nextIntervals[nextTree] - 1) > 0) {
                coalEvent[index] = nextTree;
            } else {
                coalEvent[index] = -1;
            }


            intervalLengths[index] = length;

            currentTime = nextTime;
            currentTree = nextTree;

            nextIntervals[currentTree] = nextIntervals[currentTree] + 1;

            if (nextIntervals[currentTree] == intervalsList.get(currentTree).getIntervalCount()) {
                treeFinished[currentTree] = 1;
                numTreesFinished++;
            }

            index++;  //when the loop stops, index will have value equal to the total number of intervals


        }

         System.err.println("End of while-loop");


        double[] sufficientStatistics = new double[index]; // should use smaller array.  count number of coalescent events

        int k = 0;
        double weight = 0;


        //here i indexes intervals, j indexes trees

        for (int i = 0; i < index; i++) {
            for (int j = 0; j < numTrees; j++) {
                weight = weight + lineageCounts[j][i] * (lineageCounts[j][i] - 1) * intervalLengths[i];
            }
            if (coalEvent[i] > -1) {
                sufficientStatistics[k] = weight / 2;
                k++;
                weight = 0;
            }

        }
    }

    private List<Tree> treeList;
    private List<TreeIntervals> intervalsList;

    public int nLoci() {
        return treeList.size();
    }

    public Tree getTree(int nt) {
        return treeList.get(nt);
    }

    public TreeIntervals getTreeIntervals(int nt) {
        return intervalsList.get(nt);
    }

    public double getPopulationFactor(int nt) {
        return 1.0;
    }

    public void storeTheState() {
        for (TreeIntervals intervals : intervalsList) {
            intervals.storeState();
        }
    }

    public void restoreTheState() {
        for (TreeIntervals intervals : intervalsList) {
            intervals.restoreState();
        }
    }
}


package dr.evomodel.coalescent;

import dr.evolution.coalescent.TreeIntervals;
import dr.evolution.tree.Tree;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.MatrixParameter;
import dr.inference.model.Model;
import dr.inference.model.Parameter;

import java.util.ArrayList;
import java.util.List;

import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.SymmTridiagMatrix;

/**
 * @author Mandev Gill
 * @author Marc A. Suchard
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


    // Must deal with this.  Right now, this doesn't correct anything.  The correct field length is established
    // below in calculateSufficientStatistics()
    protected int getCorrectFieldLength() {
        int tips = 0;
        for (Tree tree : treeList) {
           tips += tree.getExternalNodeCount();
        }
        //System.err.println(tips - treeList.size());
        return tips - treeList.size();
        // TODO Solution = # of internal nodes in all trees + (# of unique sampling times - 1)
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
        System.out.println("\t\tGill, Suchard, Drummond and Lemey (in preparation).");
	}

    public void wrapSetupIntervals() {
        // Do nothing
    }

    // Calculation of sufficient statistics
    //What we want:

    //-Matrix with number of lineages present in each tree during each interval
    //-length of each "w" interval and each "u" interval
    //-vector which has -1 if corresponding interval does not end with coalescent event, and otherwise has i, where i is the tree in which the coalescent event took place


    int numTrees;

    int N = 0;

    int[][] countsOfLineages;
    double[] intervalLengths;
    int[] coalEvent;

    int index;
    double[] sufficientStatistics;
    double[] coalescentIntervals;


    protected void setupSufficientStatistics() {

        index = 0;
        numTrees = treeList.size();

        for (int i = 0; i < numTrees; i++) {
            N = N + intervalsList.get(i).getIntervalCount();
        }

        countsOfLineages = new int[numTrees][N];
        intervalLengths = new double[N];
        coalEvent = new int[N];
                
        //vector keeping track of which intervals are next in each tree (in other words, how far we've gone with
        //respect to each tree
        //should be initialized to all zeros
        int nextIntervals[] = new int[numTrees];
        int treeFinished[] = new int[numTrees]; //vector indicating whether we have reached end of a certain tree
        int numTreesFinished = 0;

        //starting time in first tree
        int currentTree = 0;
        double currentTime = intervalsList.get(0).getIntervalTime(0);
        double propCurrent;

        //Figure out which tree to start in by finding smallest starting time among all trees
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

        //note that the next interval may have have same starting time as current interval
        nextIntervals[currentTree] = nextIntervals[currentTree] + 1;

        while (numTreesFinished < numTrees) {

            // Make sure that for each tree, the starting time of the next interval (nextIntervals[i]) is
            // strictly greater than currentTime.  Need this in case more than one tree has node at currentTime
            
            for(int i = 0; i < numTrees; i++){
                if(treeFinished[i] == 0){
                    while(intervalsList.get(i).getIntervalTime(nextIntervals[i]) <= currentTime){
                        nextIntervals[i] = nextIntervals[i] + 1;
                        if(nextIntervals[i] == intervalsList.get(i).getIntervalCount()){
                            treeFinished[i] = 1;
                            numTreesFinished++;
                            break;
                        }
                    }
                }
            }

            int nextTree = 0;
            double nextTime = 0;
            double length = 0;

            // In case tree 0 is finished, we look for the next unfinished tree
            while (nextTree < numTrees) {
                if (treeFinished[nextTree] == 0) {  //if this is the case should NOT move on to next tree
                    nextTime = intervalsList.get(nextTree).getIntervalTime(nextIntervals[nextTree]);
                    length = nextTime - currentTime;
                    break;
                } else {
                    nextTree++;
                }
            }

            double propNext;
            int startIndex = nextTree + 1;   //for efficiency, since we already looked at trees 0 through nextTree

            // Find tree with node closest (but not equal to) current node in terms of time
            for (int i = startIndex; i < numTrees; i++) {
                if (treeFinished[i] == 0) {
                    propNext = intervalsList.get(i).getIntervalTime(nextIntervals[i]);
                    if((propNext - currentTime) < length){
                        length = propNext - currentTime;
                        nextTime = propNext;
                        nextTree = i;
                    }
                }
            }

            // update countsOfLineages
            for (int i = 0; i < numTrees; i++) {
                if (treeFinished[i] == 0) {
                    // in case we haven't "gotten" to tree i
                    if (intervalsList.get(i).getIntervalTime(0) >= nextTime){
                        countsOfLineages[i][index] = 0;
                    }else{
                        countsOfLineages[i][index] = intervalsList.get(i).getLineageCount(nextIntervals[i]);
                        //System.err.println(i);
                        //System.err.println(index);
                        //System.err.println(nextIntervals[i]);
                        //System.err.println(countsOfLineages[i][index]);

                    }    
                }
            }

            if (intervalsList.get(nextTree).getCoalescentEvents(nextIntervals[nextTree]) > 0) {
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

        calculateSufficientStatistics();

    }




    public void calculateSufficientStatistics() {

        sufficientStatistics = new double[fieldLength];
        coalescentIntervals = new double[fieldLength];
        int m = 0;

        double weight = 0;
        double tempLength = 0;

        //here i indexes intervals, j indexes trees

        for (int i = 0; i < index; i++) {
            for (int j = 0; j < numTrees; j++) {
                weight = weight + countsOfLineages[j][i] * (countsOfLineages[j][i] - 1) * intervalLengths[i];
                }
            tempLength = tempLength + intervalLengths[i];
            if (coalEvent[i] > -1) {
                sufficientStatistics[m] = weight / 2;
                System.err.println(sufficientStatistics[m]);
                coalescentIntervals[m] = tempLength;
                m++;
                weight = 0;
                tempLength = 0;
            }

        }
    }



    
    // Likelihood implementation
    
    public double getLogLikelihood(){
        if(!likelihoodKnown){
            logCoalescentLikelihood = calculateLogCoalescentLikelihood();
            logFieldLikelihood = calculateLogFieldLikelihood();
            likelihoodKnown = true;
        }
        //System.err.println(logCoalescentLikelihood);
        //System.err.println(logFieldLikelihood);
        return logCoalescentLikelihood + logFieldLikelihood;
    }

    protected double peakLogCoalescentLikelihood() {
        return logCoalescentLikelihood;
    }

    protected double peakLogFieldLikelihood() {
        return logFieldLikelihood;
    }


    protected void setupGMRFWeights(){
        
        setupSufficientStatistics();

        //Set up the weight matrix
	    double[] offdiag = new double[fieldLength - 1];
	    double[] diag = new double[fieldLength];

	    //First set up the offdiagonal entries;

	    double rootHeight = 0;
        double prospRootHeight;

        // Set rootHeight to greatest rootHeight among all trees.  Make sure this is what we want.
        for(int i = 0; i < numTrees; i++){
            prospRootHeight = treeList.get(i).getNodeHeight(treeList.get(i).getRoot());
            if(prospRootHeight > rootHeight){
                rootHeight = prospRootHeight;
            }
        }


	    for (int i = 0; i < (fieldLength - 1); i++) {
		    offdiag[i] = (-2.0 / (coalescentIntervals[i] + coalescentIntervals[i + 1])) * rootHeight;
	    }


	    //Then set up the diagonal entries;
	    for (int i = 1; i < (fieldLength - 1); i++)
	    diag[i] = -(offdiag[i] + offdiag[i - 1]);

	    //Take care of the endpoints
	    diag[0] = -offdiag[0];
	    diag[fieldLength - 1] = -offdiag[fieldLength - 2];

	    weightMatrix = new SymmTridiagMatrix(diag, offdiag);
    }




    public SymmTridiagMatrix getScaledWeightMatrix(double precision) {
		SymmTridiagMatrix a = weightMatrix.copy();
		for (int i = 0; i < a.numRows() - 1; i++) {
			a.set(i, i, a.get(i, i) * precision);
			a.set(i + 1, i, a.get(i + 1, i) * precision);
		}
		a.set(fieldLength - 1, fieldLength - 1, a.get(fieldLength - 1, fieldLength - 1) * precision);
		return a;
	}


    public SymmTridiagMatrix getStoredScaledWeightMatrix(double precision) {
		SymmTridiagMatrix a = storedWeightMatrix.copy();
		for (int i = 0; i < a.numRows() - 1; i++) {
			a.set(i, i, a.get(i, i) * precision);
			a.set(i + 1, i, a.get(i + 1, i) * precision);
		}
		a.set(fieldLength - 1, fieldLength - 1, a.get(fieldLength - 1, fieldLength - 1) * precision);
		return a;
	}

	public SymmTridiagMatrix getScaledWeightMatrix(double precision, double lambda) {
		if (lambda == 1)
			return getScaledWeightMatrix(precision);

		SymmTridiagMatrix a = weightMatrix.copy();
		for (int i = 0; i < a.numRows() - 1; i++) {
			a.set(i, i, precision * (1 - lambda + lambda * a.get(i, i)));
			a.set(i + 1, i, a.get(i + 1, i) * precision * lambda);
		}

		a.set(fieldLength - 1, fieldLength - 1, precision * (1 - lambda + lambda * a.get(fieldLength - 1, fieldLength - 1)));
		return a;
	}


    /*
	public double[] getCoalescentIntervalHeights() {
		double[] a = new double[coalescentIntervals.length];

		a[0] = coalescentIntervals[0];

		for (int i = 1; i < a.length; i++) {
			a[i] = a[i - 1] + coalescentIntervals[i];
		}
		return a;
	}

	public SymmTridiagMatrix getCopyWeightMatrix() {
		return weightMatrix.copy();
	}
	*/

	public SymmTridiagMatrix getStoredScaledWeightMatrix(double precision, double lambda) {
		if (lambda == 1)
			return getStoredScaledWeightMatrix(precision);

		SymmTridiagMatrix a = storedWeightMatrix.copy();
		for (int i = 0; i < a.numRows() - 1; i++) {
			a.set(i, i, precision * (1 - lambda + lambda * a.get(i, i)));
			a.set(i + 1, i, a.get(i + 1, i) * precision * lambda);
		}

		a.set(fieldLength - 1, fieldLength - 1, precision * (1 - lambda + lambda * a.get(fieldLength - 1, fieldLength - 1)));
		return a;
	}

    /*
	protected void storeState() {
		super.storeState();
		System.arraycopy(coalescentIntervals, 0, storedCoalescentIntervals, 0, coalescentIntervals.length);
		System.arraycopy(sufficientStatistics, 0, storedSufficientStatistics, 0, sufficientStatistics.length);
		storedWeightMatrix = weightMatrix.copy();
        storedLogFieldLikelihood = logFieldLikelihood;
	}


	protected void restoreState() {
		super.restoreState();
		System.arraycopy(storedCoalescentIntervals, 0, coalescentIntervals, 0, storedCoalescentIntervals.length);
		System.arraycopy(storedSufficientStatistics, 0, sufficientStatistics, 0, storedSufficientStatistics.length);
		weightMatrix = storedWeightMatrix;
        logFieldLikelihood = storedLogFieldLikelihood;

	}

	protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type){
		likelihoodKnown = false;
        // Parameters (precision and popsizes do not change intervals or GMRF Q matrix
	}
    */
    


    private double calculateLogCoalescentLikelihood() {

		if (!intervalsKnown) {
			// intervalsKnown -> false when handleModelChanged event occurs in super.
			wrapSetupIntervals();
			setupGMRFWeights();
		}

		// Matrix operations taken from block update sampler to calculate data likelihood and field prior

		double currentLike = 0;
        double[] currentGamma = popSizeParameter.getParameterValues();

		for (int i = 0; i < fieldLength; i++) {
			currentLike += -currentGamma[i] - sufficientStatistics[i] * Math.exp(-currentGamma[i]);
		}

		return currentLike;// + LogNormalDistribution.logPdf(Math.exp(popSizeParameter.getParameterValue(coalescentIntervals.length - 1)), mu, sigma);
	}

    private double calculateLogFieldLikelihood() {

        if (!intervalsKnown) {
            // intervalsKnown -> false when handleModelChanged event occurs in super.
            wrapSetupIntervals();
            setupGMRFWeights();
        }

        double currentLike = 0;
        DenseVector diagonal1 = new DenseVector(fieldLength);
        DenseVector currentGamma = new DenseVector(popSizeParameter.getParameterValues());

        SymmTridiagMatrix currentQ = getScaledWeightMatrix(precisionParameter.getParameterValue(0), lambdaParameter.getParameterValue(0));
        currentQ.mult(currentGamma, diagonal1);

        currentLike += 0.5 * (fieldLength - 1) * Math.log(precisionParameter.getParameterValue(0)) - 0.5 * currentGamma.dot(diagonal1);
        if (lambdaParameter.getParameterValue(0) == 1) {
            currentLike -= (fieldLength - 1) / 2.0 * LOG_TWO_TIMES_PI;
        } else {
            currentLike -= fieldLength / 2.0 * LOG_TWO_TIMES_PI;
        }

        return currentLike;
    }



    private double logCoalescentLikelihood;
    private double logFieldLikelihood;
    private boolean likelihoodKnown = false;



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


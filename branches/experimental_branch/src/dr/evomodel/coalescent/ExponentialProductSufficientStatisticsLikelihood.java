package dr.evomodel.coalescent;

import java.util.List;

import dr.evolution.coalescent.IntervalType;
import dr.evolution.tree.Tree;
import dr.evomodel.coalescent.OldAbstractCoalescentLikelihood.CoalescentEventType;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Likelihood;

/**
 * Calculates a product of exponential densities and exponential tail probabilities.
 *
 * @author Guy Baele
 */

public class ExponentialProductSufficientStatisticsLikelihood extends OldAbstractCoalescentLikelihood {
	
	//not used at the moment
	public final static boolean FIXED_TREE = false;
	
	private TreeModel treeModel;
	private double[] posteriorMeans;
	
	protected int fieldLength;
	protected double[] coalescentIntervals;
    protected double[] sufficientStatistics;

	//make sure to use in combination with coalescentEventsStatistic
	public ExponentialProductSufficientStatisticsLikelihood(TreeModel treeModel, double[] posteriorMeans) {
		super("ExponentialProductSufficientStatisticsLikelihood");
		this.treeModel = treeModel;
		this.posteriorMeans = posteriorMeans;
		
		tree = treeModel;
		addModel((TreeModel) tree);
		
		this.fieldLength = posteriorMeans.length;
		wrapSetupIntervals();
		coalescentIntervals = new double[fieldLength];
        sufficientStatistics = new double[fieldLength];
	}
	
	protected void wrapSetupIntervals() {
        setupIntervals();
    }
	
	/*protected void setTree(List<Tree> treeList) {
        if (treeList.size() != 1) {
            throw new RuntimeException("GMRFSkyrideLikelihood only implemented for one tree");
        }
        this.tree = treeList.get(0);
        this.treesSet = null;
        if (tree instanceof TreeModel) {
            addModel((TreeModel) tree);
        }
    }*/
	
	protected void setupSufficientStatistics() {
        int index = 0;

        double length = 0;
        double weight = 0;
        for (int i = 0; i < getIntervalCount(); i++) {
        	//System.err.println(getInterval(i));
            length += getInterval(i);
            weight += getInterval(i) * getLineageCount(i) * (getLineageCount(i) - 1);
            if (getIntervalType(i) == CoalescentEventType.COALESCENT) {
                coalescentIntervals[index] = length;
                sufficientStatistics[index] = weight / 2.0;
                index++;
                length = 0;
                weight = 0;
            }
        }
    }
	
	private void makeIntervalsKnown() {
        if (!intervalsKnown) {
            wrapSetupIntervals();
            intervalsKnown = true;
        }
    }
	
	public double calculateLogLikelihood() {
		tree = treeModel;
        makeIntervalsKnown();
        setupSufficientStatistics();

        // Matrix operations taken from block update sampler to calculate data likelihood and field prior

        double currentLike = 0;
        double[] currentGamma = posteriorMeans;

        for (int i = 0; i < fieldLength; i++) {
        	//System.err.println(currentGamma[i] + "    " + sufficientStatistics[i] + "    " + (-currentGamma[i] - sufficientStatistics[i] * Math.exp(-currentGamma[i])));
            currentLike += -currentGamma[i] - sufficientStatistics[i] * Math.exp(-currentGamma[i]);
        }
        //System.err.println("currentLike = " + currentLike + "\n");

        return currentLike;
    }

	/**
	 * Overridden to always return false.
	 */
	protected boolean getLikelihoodKnown() {
		return false;
	}

}

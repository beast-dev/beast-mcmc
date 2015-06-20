package dr.evomodel.coalescent;

import dr.evomodel.tree.TreeModel;
import dr.inference.model.Likelihood;

/**
 * Calculates a product of exponential densities and exponential tail probabilities.
 *
 * @author Guy Baele
 */

public class ExponentialProductPosteriorMeansLikelihood extends Likelihood.Abstract {
	
	//not used at the moment
	public final static boolean FIXED_TREE = false;
    public static final boolean DEBUG = true;
	
	private TreeModel treeModel;
	private double[] posteriorMeans;

	//make sure to use in combination with coalescentEventsStatistic
	public ExponentialProductPosteriorMeansLikelihood(TreeModel treeModel, double[] posteriorMeans) {
		super(treeModel);
		this.treeModel = treeModel;
		this.posteriorMeans = posteriorMeans;
	}
	
	public double calculateLogLikelihood() {
		
		//System.err.println(treeModel);
		
		double logPDF = 0.0;
		
		//System.err.println("log(popSize) = " + this.popSize);
		
		CoalescentTreeIntervalStatistic ctis = new CoalescentTreeIntervalStatistic(treeModel);
		int coalescentIntervalCounter = 0;

        if (DEBUG) {
            System.err.println("ExponentialProductPosteriorMeansLikelihood dimension: " + ctis.getDimension());
        }
		
		for (int i = 0; i < ctis.getDimension(); i++) {
			
			int combinations = (int)ctis.getLineageCount(i)*((int)ctis.getLineageCount(i)-1)/2;
			double branchLength = ctis.getStatisticValue(i);
			
			//System.err.println("combinations = " + combinations);
			//System.err.println("branchLength = " + branchLength);
			
			//System.err.println(ctis.getLineageCount(i));
			
			//single-lineage intervals are not counted
			if (ctis.getLineageCount(i) != 1) {
				
				//System.err.println(i + " -> lineage count: " + ctis.getLineageCount(i));
				
				if (i == ctis.getDimension()-1) {
					//coalescent event at root: exponential density
					//System.err.print("coalescent event at root: ");
					double logContribution = -posteriorMeans[coalescentIntervalCounter] - combinations*branchLength*Math.exp(-posteriorMeans[coalescentIntervalCounter]);
					if (DEBUG) {
                        System.err.println(i + ": " + logContribution);
                    }
                    logPDF += logContribution;
					//System.err.println(logContribution);
				} else if (ctis.getLineageCount(i) > ctis.getLineageCount(i+1)) {
					//coalescent event: exponential density
					//System.err.print("coalescent event (not at root): ");
					double logContribution = -posteriorMeans[coalescentIntervalCounter] - combinations*branchLength*Math.exp(-posteriorMeans[coalescentIntervalCounter]);
                    if (DEBUG) {
                        System.err.println(i + ": " + logContribution);
                    }
                    logPDF += logContribution;
					//System.err.println(logContribution);
					coalescentIntervalCounter++;
				} else {
					//sampling event: exponential tail probability
					//System.err.print("sampling event: ");
					double logContribution = -combinations*branchLength*Math.exp(-posteriorMeans[coalescentIntervalCounter]);
                    if (DEBUG) {
                        System.err.println(i + ": " + logContribution);
                    }
                    logPDF += logContribution;
					//System.err.println(logContribution);
				}
				
			}
			
		}
		
		//System.err.println("expoLike = " + logPDF + "\n");
		
		return logPDF;
		
	}

	/**
	 * Overridden to always return false.
	 */
	protected boolean getLikelihoodKnown() {
		return false;
	}

}

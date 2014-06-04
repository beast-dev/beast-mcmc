package dr.evomodel.coalescent;

import dr.evomodel.tree.TreeModel;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.math.distributions.GammaDistribution;

/**
 * Calculates a product of exponential densities and exponential tail probabilities.
 *
 * @author Guy Baele
 */

public class ExponentialProductLikelihood extends Likelihood.Abstract {
	
	public final static boolean REDUCE_TO_EXPONENTIAL = true;
	
	//not used at the moment
	public final static boolean FIXED_TREE = false;
	
	private TreeModel treeModel;
	private double popSize;

	//make sure to provide a log(popSize)
	public ExponentialProductLikelihood(TreeModel treeModel, double popSize) {
		super(treeModel);
		this.treeModel = treeModel;
		this.popSize = popSize;
	}
	
	public double calculateLogLikelihood() {
		
		//System.err.println(treeModel);
		
		double logPDF = 0.0;
		
		//System.err.println("log(popSize) = " + this.popSize);
		
		CoalescentTreeIntervalStatistic ctis = new CoalescentTreeIntervalStatistic(treeModel);
		
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
					double logContribution = -popSize - combinations*branchLength*Math.exp(-popSize);
					logPDF += logContribution;
					//System.err.println(logContribution);
				} else if (ctis.getLineageCount(i) > ctis.getLineageCount(i+1)) {
					//coalescent event: exponential density
					//System.err.print("coalescent event (not at root): ");
					double logContribution = -popSize - combinations*branchLength*Math.exp(-popSize);
					logPDF += logContribution;
					//System.err.println(logContribution);
				} else {
					//sampling event: exponential tail probability
					//System.err.print("sampling event: ");
					double logContribution = -combinations*branchLength*Math.exp(-popSize);
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

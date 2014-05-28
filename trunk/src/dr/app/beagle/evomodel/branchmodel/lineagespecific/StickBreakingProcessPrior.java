package dr.app.beagle.evomodel.branchmodel.lineagespecific;

import dr.app.beagle.evomodel.branchmodel.BranchModel;
import dr.app.bss.Utils;
import dr.evolution.tree.NodeRef;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.math.GammaFunction;
import dr.math.distributions.Distribution;
import dr.math.distributions.MultivariateDistribution;

public class StickBreakingProcessPrior implements MultivariateDistribution {

	private Distribution baseDistribution;
	private Parameter mdParameter;
	private int categoryCount;
	
	private TreeModel treeModel;
	private BranchModel branchSpecificModel;
	
	public StickBreakingProcessPrior(
			BranchModel branchSpecificModel, //
			TreeModel treeModel, //
			Distribution baseDistribution,//
//			List<Parameter> parameters
			Parameter mdParameter
			) {
		
		this.baseDistribution = baseDistribution;
		this.mdParameter = mdParameter;
		
		this.branchSpecificModel = branchSpecificModel;
		this.treeModel = treeModel;
		
		this.categoryCount = mdParameter.getDimension();
		
	}//END: Constructor
	
	
	private int[] getCounts() {

		int[] counts = new int[categoryCount];
		
		for (NodeRef branch : treeModel.getNodes()) {

			int branchParameterIndex = ((BranchSpecific) branchSpecificModel)
					.getBranchModelMapping(branch).getOrder()[0];
			counts[branchParameterIndex]++;
		}//END: branch loop

		return counts;
	}// END: getBranchAssignmentCounts
	
	
	private double[] getValues() {
		
		double[] values = new double[categoryCount];
		

		for(int i=0;i<categoryCount;i++) {
			
			values[i]=mdParameter.getParameterValue(i);
			
		}
		
		
		return values;
	}//END: getValues
	
	/*
	 * Distribution Likelihood
	 */
	public double getLogLikelihood() {

		double logLike = 0.0;

		int[] counts = getCounts();
		double[] values = getValues();
		double countSum = Utils.sumArray(counts);
		
		logLike += GammaFunction.lnGamma(countSum);
		for (int i = 0; i < categoryCount; i++) {

			logLike += ( (counts[i] - 1) * Math.log(values[i]) - GammaFunction.lnGamma(counts[i])  ); 
					
		}// END: i loop

		return logLike;
	}// END: getLogLikelihood


	@Override
	public double logPdf(double[] x) {
		double logLike = 0.0;

		int[] counts = getCounts();
		double countSum = Utils.sumArray(counts);
		
		logLike += GammaFunction.lnGamma(countSum);
		for (int i = 0; i < categoryCount; i++) {

			logLike += (  (counts[i] - 1) * Math.log(x[i]) -  GammaFunction.lnGamma(counts[i])  ); 
					
		}// END: i loop

		return logLike;
	}//END: logPdf

	@Override
	public double[] getMean() {
        
		int[] counts = getCounts();
		double countSum = Utils.sumArray(counts);
		
		double[] mean = new double[categoryCount];
        
        for (int i = 0; i < categoryCount; i++) {
            mean[i] = counts[i] / countSum;
        }
        
        return mean;
	}//END: mean

	@Override
	public double[][] getScaleMatrix() {
		return null;
	}

	@Override
	public String getType() {
		return null;
	}
	
}// END: class


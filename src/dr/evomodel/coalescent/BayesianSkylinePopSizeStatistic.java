package dr.evomodel.coalescent;

import dr.inference.model.Statistic;

public class BayesianSkylinePopSizeStatistic extends Statistic.Abstract{

	public double time;
	public BayesianSkylineLikelihood bsl;

	public BayesianSkylinePopSizeStatistic(double time, 
			BayesianSkylineLikelihood bsl){
		this.time = time;
		this.bsl = bsl;
	}
	
	public int getDimension() {
		return 1;
	}

	public double getStatisticValue(int dim) {
		double[] heights = bsl.getGroupHeights();
		double[] sizes = bsl.getPopSizeParameter().getParameterValues();
		
		
		for(int i = 0; i < heights.length; i++){
			if(this.time < heights[i]){
				return sizes[i];
			}
		}
		
		return Double.NaN;
	}

}

package dr.evomodel.coalescent;

import dr.inference.model.Statistic;

public class GMRFPopSizeStatistic extends Statistic.Abstract{
	
	private GMRFSkyrideLikelihood gsl;
	private double[] time;

	public GMRFPopSizeStatistic(double[] time, GMRFSkyrideLikelihood gsl){
		super("Popsize");
		this.gsl = gsl;
		this.time = time;
	}
	
	public String getDimensionName(int i){
		return getStatisticName() + Double.toString(time[i]);
	}
	
	public int getDimension() {
		return time.length;
	}

	public double getStatisticValue(int dim) {
		double[] coalescentHeights = gsl.getCoalescentIntervalHeights();
		double[] popSizes = gsl.getPopSizeParameter().getParameterValues();
		
		assert popSizes.length == coalescentHeights.length;
		
		for(int i = 0; i < coalescentHeights.length; i++){
			if(time[dim] < coalescentHeights[i]){
				return popSizes[i];
			}
		}
				
		return Double.NaN;
	}


}

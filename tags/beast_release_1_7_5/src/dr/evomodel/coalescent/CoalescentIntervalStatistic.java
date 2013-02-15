package dr.evomodel.coalescent;

import dr.inference.model.Statistic;


public class CoalescentIntervalStatistic extends Statistic.Abstract{

	private GMRFSkyrideLikelihood acl;
	private int dimension;
	
	public CoalescentIntervalStatistic(GMRFSkyrideLikelihood acl){
		this.acl = acl;
		
		dimension = acl.getCoalescentIntervalHeights().length;
	}
	
	public int getDimension() {
		return dimension;
	}

	public double getStatisticValue(int dim) {
			return acl.getCoalescentIntervalHeights()[dim];
	}

}

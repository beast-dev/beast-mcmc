package dr.evomodel.coalescent;

import dr.inference.model.Statistic;


public class CoalescentIntervalStatistic extends Statistic.Abstract {
	
	private coalescentInterval coalInt;
	
	public CoalescentIntervalStatistic(GMRFSkyrideLikelihood acl) {
		this.coalInt = new GMRFStatistic(acl);
	}
	
	public CoalescentIntervalStatistic(CoalescentLikelihood coal) {
		this.coalInt = new CoalescentStatistic(coal);
	}
	
	public int getDimension() {
		return coalInt.getDimension();
	}

	public double getStatisticValue(int dim) {
		return coalInt.getStatisticValue(dim);
	}
	
	private class GMRFStatistic implements coalescentInterval {
		
		private GMRFSkyrideLikelihood acl;
		private int dimension;
		
		private GMRFStatistic(GMRFSkyrideLikelihood acl) {
			this.acl = acl;
			this.dimension = acl.getCoalescentIntervalHeights().length;
		}
		
		public int getDimension() {
			return dimension;
		}
		
		public double getStatisticValue(int dim) {
			return acl.getCoalescentIntervalHeights()[dim];
		}
		
	}
	
	private class CoalescentStatistic implements coalescentInterval {
		
		private CoalescentLikelihood coal;
		private int dimension;
		
		private CoalescentStatistic(CoalescentLikelihood coal) {
			this.coal = coal;
			this.dimension = coal.getCoalescentIntervalHeights().length;
		}
		
		public int getDimension() {
			return dimension;
		}
		
		public double getStatisticValue(int dim) {
			return coal.getCoalescentIntervalHeights()[dim];
		}
		
	}
	
	private interface coalescentInterval {
		
		public int getDimension();
		
		public double getStatisticValue(int dim);
		
	}

}

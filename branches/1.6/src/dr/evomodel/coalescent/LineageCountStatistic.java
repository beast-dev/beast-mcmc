package dr.evomodel.coalescent;

import dr.inference.model.Statistic;

public class LineageCountStatistic extends Statistic.Abstract{
	private OldAbstractCoalescentLikelihood acl;
	
	public LineageCountStatistic(OldAbstractCoalescentLikelihood acl){
		this.acl = acl;
	}
	
	public int getDimension() {
		return acl.getIntervalCount();
	}

	public double getStatisticValue(int dim) {
			return acl.getLineageCount(dim);
	}

}

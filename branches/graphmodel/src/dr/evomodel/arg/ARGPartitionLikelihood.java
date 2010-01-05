package dr.evomodel.arg;

import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.Model;

public abstract class ARGPartitionLikelihood extends AbstractModelLikelihood {

	
	private final ARGModel arg;
	
	public abstract double[] generatePartition();
	public abstract double getLogLikelihood(double[] partition);
	
	public double getLogLikelihood() {
    			
    	double logPrior = 0;
    	
    	for(int i = 0, n = getReassortmentNodeCount(); i < n; i++){
    		logPrior += getLogLikelihood(getPartition(i));
    	}
    	
    	return logPrior;
	}
	
	public double[] getPartition(int i){
		if(arg.getReassortmentNodeCount() == 0){
			return null;
		}
		
		return arg.getPartitioningParameters().getParameter(i).getParameterValues();
	}
	
	public ARGPartitionLikelihood(String id, ARGModel arg){
		super(id);
		
		
		this.arg = arg;
	}
	
	public int getNumberOfPartitionsMinusOne(){
		return arg.getNumberOfPartitions() - 1;
	}
	
	public int getReassortmentNodeCount(){
		return arg.getReassortmentNodeCount();
	}
	
	public Model getModel() {
		return this;
	}
	
	public void makeDirty() {
		
	}
}

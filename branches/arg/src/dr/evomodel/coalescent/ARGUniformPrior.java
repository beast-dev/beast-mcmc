package dr.evomodel.coalescent;

import dr.evomodel.tree.ARGModel;

public class ARGUniformPrior extends CoalescentLikelihood{

	public static final String ARG_UNIFORM_PRIOR = "argUniformPrior";
	
	private ARGModel arg;
	
	ARGUniformPrior(String name, ARGModel arg) {
		super(ARG_UNIFORM_PRIOR);
		
		this.arg = arg;
	}
	
	public double getLogLikelihood(){
		if(likelihoodKnown)
			return logLikelihood;
		
		likelihoodKnown = true;
		logLikelihood = calculateLogLikelihood();
		
		return logLikelihood;
	}
	
	public double calculateLogLikelihood(){
		
		double treeHeight = arg.getNodeHeight(arg.getRoot());
		int internalNodes = arg.getInternalNodeCount() - 1;
		
		double logLike = logFactorial(internalNodes) - (double)internalNodes*Math.log(treeHeight);
		
		if(arg.getReassortmentNodeCount() == 0){
			logLike -= Math.log(18);
		}else{
			
		}
		
		return 0.0;
	}
	
	public void storeState(){
		 storedLikelihoodKnown = likelihoodKnown;
	     storedLogLikelihood = logLikelihood;
	}
	
	public void restoreState(){
		likelihoodKnown = storedLikelihoodKnown;
        logLikelihood = storedLogLikelihood;
	}
	
	private double logFactorial(int n){
		double rValue = 0;
		
		for(int i = n; i > 0; i--){
			rValue += Math.log(rValue);
		}
		return rValue; 
	}
	

	
	
}

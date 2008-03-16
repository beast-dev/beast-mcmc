package dr.evomodel.operators;

import java.util.logging.Logger;

import dr.evomodel.tree.ARGModel;
import dr.inference.model.VariableSizeCompoundParameter;
import dr.inference.operators.OperatorFailedException;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;

public class ARGReassortmentOperator extends SimpleMCMCOperator{

	public static final String ADD_PROBABILITY = "addProbability";
	public static final String ARG_REASSORTMENT_OPERATOR = "argReassortmentOperator";
	public static final String INTERNAL_NODES = "internalNodes";
	public static final String INTERNAL_AND_ROOT = "internalNodesPlusRoot";
	public static final String NODE_RATES = "nodeRates";
	public static final String CHOOSE_BRANCHES_FIRST = "chooseBranchesFirst";
	public static final double LOG_TWO = Math.log(2.0);
	
	private double singlePartitionProbability = 0.0;
	private double probBelowRoot = 0.9;
	private double size = 0.0; //Used to choose add step
	private boolean branchesFirst;
	
	private ARGModel arg;
	
	private VariableSizeCompoundParameter internalNodeParameters;
	private VariableSizeCompoundParameter internalAndRootNodeParameters;
	private VariableSizeCompoundParameter nodeRates;
	
	public ARGReassortmentOperator(ARGModel arg, int weight, boolean branchesFirst,
								   double singlePartProb, double addProbability,
								   VariableSizeCompoundParameter internalNodeParameters,
								   VariableSizeCompoundParameter internalAndRootNodeParameters,
							   	   VariableSizeCompoundParameter nodeRates){
		this.arg = arg;
		this.internalNodeParameters = internalNodeParameters;
		this.internalAndRootNodeParameters = internalAndRootNodeParameters;
		this.nodeRates = nodeRates;
		
		this.branchesFirst = branchesFirst;
		this.singlePartitionProbability = singlePartProb;
		
		setWeight(weight);
		
		this.size = addProbability;
		
		//Transformed for computation reasons
		probBelowRoot = -Math.log(1 - Math.sqrt(probBelowRoot)); 
	}
	
	
	public double doOperation() throws OperatorFailedException {
		double logHastings = 0;
		
		try{
			if (MathUtils.nextDouble() < 1.0/(1 + Math.exp(-size)))
				logHastings = addOperation() - size;
			else
				logHastings = removeOperation() + size;
		}catch(NoReassortmentEventException nree){
			return Double.NEGATIVE_INFINITY;
		}catch(OperatorFailedException ofe){
			Logger.getLogger("dr.evomodel.operators").fine(ofe.getMessage());
		}
		
		return logHastings;
	}

	private double addOperation() throws OperatorFailedException{
		if(branchesFirst)
			return addOperationBranchesFirst();
					
		return addOperationHeightsFirst();
	}
	
	private double addOperationBranchesFirst() throws OperatorFailedException{
		
		double logHastings = 0;
		double treeHeight = arg.getNodeHeight(arg.getRoot());
		double newBifurcationHeight  = Double.POSITIVE_INFINITY; 
		double newReassortmentHeight = Double.POSITIVE_INFINITY;
		
		double theta = probBelowRoot / treeHeight;
		
		while(newBifurcationHeight > treeHeight && newReassortmentHeight > treeHeight){
			newBifurcationHeight = MathUtils.nextExponential(theta);
			newReassortmentHeight = MathUtils.nextExponential(theta);
		}
		
		logHastings += theta*(newBifurcationHeight + newReassortmentHeight) - LOG_TWO
		 			- 2.0*Math.log(theta) + Math.log(1 - Math.exp(-2.0*treeHeight*theta));
		
		if(newBifurcationHeight < newReassortmentHeight){
			double temp = newBifurcationHeight;
			newBifurcationHeight = newReassortmentHeight;
			newReassortmentHeight = temp;
		}
		
		return 0;
	}
	
	private double addOperationHeightsFirst() throws OperatorFailedException{
		
		return 0;
	}
	
	private double removeOperation() throws OperatorFailedException{
	
		return 0;
	}
	
	public String getOperatorName() {
		return ARG_REASSORTMENT_OPERATOR;
	}

	public String getPerformanceSuggestion() {
		return "Try changing the add probability probability";
	}

	private class NoReassortmentEventException extends OperatorFailedException {
		private static final long serialVersionUID = 1L;
		public NoReassortmentEventException() {
			super("");
		}
	}
}

package dr.inference.operators;


	import java.util.List;

import dr.inference.model.Parameter;
	import dr.math.MathUtils;
import dr.inference.model.VectorSliceParameter;
	/**
	 * An operator that flips bits in same position for all strata and hieararchical level.
	 *
	 * 
	 *  @author Gabriela Cybis
	 */
	public class HierarchicalBitFlipOperator extends SimpleMCMCOperator {

	    public HierarchicalBitFlipOperator(Parameter hParameter, Parameter[] strataParameters, int NStrata, double weight, boolean usesPriorOnSum) {

	    	this.Nstrata = NStrata;
	    	this.hParameter = hParameter;
	    	this.strataParameters = strataParameters;
	    	this.usesPriorOnSum = usesPriorOnSum;
	        setWeight(weight);
	        
	    }

	   
	    

	    /**
	     * Change the parameter and return the hastings ratio.
	     * Flip (Switch a 0 to 1 or 1 to 0) for a random bit in a bit vector.
	     * Return the hastings ratio which makes all subsets of vectors with the same number of 1 bits
	     * equiprobable
	     */
	    public final double doOperation() {
	        final int dim = hParameter.getDimension();
	        double logq = 0.0;
	       
	        //comeca novo
	        int rep= MathUtils.nextInt(4);
	        
	        for (int repete=0; repete<rep+1;repete++ ){
	    //    double sum = 0.0;

	        if(usesPriorOnSum) {
	            for (int i = 0; i < dim; i++) {
	      //          sum += hParameter.getParameterValue(i);
	            }
	        }

	      
	        
	        final int pos = MathUtils.nextInt(dim);

	        int value = (int) hParameter.getParameterValue(pos);
	        
	        if (value == 0) {
	            hParameter.setParameterValue(pos, 1.0);

	        //    if(usesPriorOnSum)
	          //      logq = -Math.log((dim - sum) / (sum + 1));

	        } else if (value == 1) {
	            hParameter.setParameterValue(pos, 0.0);
	         //   if(usesPriorOnSum)
	           //     logq = -Math.log(sum / (dim - sum + 1));

	        } else {
	            throw new RuntimeException("expected 1 or 0");
	        }

	        for (int j=0; j<Nstrata ; j++){
	        	
	        	value = (int) strataParameters[j].getParameterValue(pos);
	 	       
	 	        if (value == 0) {
	 	            strataParameters[j].setParameterValue(pos, 1.0);

	 	        } else if (value == 1) {
	 	            strataParameters[j].setParameterValue(pos, 0.0);
	 	            
	 	        } else {
	 	            throw new RuntimeException("expected 1 or 0");
	 	        }
	 
	        
	        }
	        }  
	    
	        return logq;
	    }

	    // Interface MCMCOperator
	    public final String getOperatorName() {
	        return "bitFlip(" + hParameter.getParameterName() + ")";
	    }

	    public final String getPerformanceSuggestion() {
	        return "no performance suggestion";
	    }

	    public String toString() {
	        return getOperatorName();
	    }

	    // Private instance variables

	    private Parameter hParameter = null;
	    private boolean usesPriorOnSum = false;
	    private Parameter[] strataParameters = null ;
	    int Nstrata; 
	}



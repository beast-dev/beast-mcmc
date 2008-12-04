package dr.evomodel.arg;

import dr.evomodel.arg.operators.ARGPartitioningOperator;
import dr.inference.model.Likelihood;
import dr.math.MathUtils;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLObjectParser;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;


public class UniformPartitionLikelihood extends ARGPartitionLikelihood implements Likelihood {

    public static final String UNIFORM_PARTITION_LIKELIHOOD = "uniformPartitionLikelihood";
    private double logStateCount;
    private boolean isRecombination;
	
    
    public UniformPartitionLikelihood(String id, ARGModel arg){
    	super(id,arg);
    	
    	isRecombination = arg.isRecombinationPartitionType();
    	
    	int numberOfPartitionsMinusOne = getNumberOfPartitionsMinusOne();
    	
    	if(arg.isRecombinationPartitionType()){
    		//For example, if we have five partitions, we can have the following
    		// 0 1 1 1 1
    		// 0 0 1 1 1       <-- Equates to four possibilities
    		// 0 0 0 1 1   
    		// 0 0 0 0 1
    		
    		logStateCount = Math.log(numberOfPartitionsMinusOne);
    	}else{
    		//You basically choose a subset from all possible subset of the final four
        	//there are 2^numberOfPartitionsMinusOne of these
        	//
        	//Except! you cannot choose the empty set.  
        	
    		
        	double rValue = Math.pow(2.0,numberOfPartitionsMinusOne) - 1;
        	
        	logStateCount = Math.log(rValue);
    	}
    }
    
    public double[] generatePartition(){
    	if(isRecombination){
    		return generateRecombinationPartition();
    	}
    	return generateReassortmentPartition();
    }
    
    private double[] generateRecombinationPartition(){
    	int numberOfPartitionsMinusOne = getNumberOfPartitionsMinusOne();

        int cut = MathUtils.nextInt(numberOfPartitionsMinusOne);

        int leftValue = 0;  //At one time, these values could switch.
        int rightValue = 1;

        double[] partition = new double[numberOfPartitionsMinusOne + 1];
        
        for (int i = 0; i < cut + 1; i++)
            partition[i] = leftValue;
        for (int i = cut + 1; i < partition.length; i++)
            partition[i] = rightValue;
        
        return partition;
    }
    
    private double arraySum(double[] x){
    	double a = 0;
    	for(double b : x)
    		a += b;
    	return a;
    }
    
    private double[] generateReassortmentPartition(){
    	int numberOfPartitions = getNumberOfPartitionsMinusOne() + 1;

        double[] partition = new double[numberOfPartitions];

        while (arraySum(partition) == 0) {
            for (int i = 1; i < partition.length; i++) {
                if (MathUtils.nextBoolean()) {
                    partition[i] = 1.0;
                } else {
                    partition[i] = 0.0;
                }
            }
        }
        
        partition[0] = 0.0;       
       
        return partition;
    }
    
    
    public double getLogLikelihood(double[] partition){
    	return -logStateCount;
    }
    
    public static XMLObjectParser PARSER = new AbstractXMLObjectParser(){

		public String getParserDescription() {
			return null;
		}

		public Class getReturnType() {
			return UniformPartitionLikelihood.class;
		}

		public XMLSyntaxRule[] getSyntaxRules() {
			return new XMLSyntaxRule[]{
				new ElementRule(ARGModel.class, false),	
			};
		}

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {
			ARGModel arg = (ARGModel)xo.getChild(ARGModel.class);
			
			String id = ""; 
			if(xo.hasId())
				id = xo.getId();
			
			return new UniformPartitionLikelihood(id,arg);
		}

		public String getParserName() {
			return UNIFORM_PARTITION_LIKELIHOOD;
		}
    	
    };

}

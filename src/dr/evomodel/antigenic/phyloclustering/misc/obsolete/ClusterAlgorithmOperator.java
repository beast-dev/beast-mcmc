package dr.evomodel.antigenic.phyloclustering.misc.obsolete;

import dr.inference.model.MatrixParameter;
import dr.inference.model.Parameter;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.distributions.MultivariateNormalDistribution;
import dr.xml.*;


/**
 * A Gibbs operator for allocation of items to clusters under a distance dependent Chinese restaurant process.
 *
 * @author Charles Cheung
 * @author Trevor Bedford
 */
public class ClusterAlgorithmOperator extends SimpleMCMCOperator  {

	//Parameter locationDrift;  // no longer need to know
	Parameter virusOffsetsParameter;
	private double sigmaSq =1;
	private int numdata = 0; //NEED TO UPDATE
	//private double[] groupSize; 
    
    private MatrixParameter mu = null;
    private Parameter clusterLabels = null;
    private Parameter K = null;
    
    private MatrixParameter virusLocations = null;
    
    private int maxLabel = 0;
    private int[] muLabels = null;

    private int[] groupSize;
  //  public ClusterViruses clusterLikelihood = null;


   private double numAcceptMoveMu = 0;
   private double numProposeMoveMu = 0;
   
   private double numAcceptMoveC = 0;
   private double numProposeMoveC = 0;
   
   private int isMoveMu = -1;
    
   
	private double[] old_vLoc0 ;
	private double[] old_vLoc1 ;

    private Parameter clusterOffsetsParameter;
   
    private int groupSelectedChange = -1;
    private int virusIndexChange = -1;
    private double originalValueChange = -1;
    private int dimSelectChange = -1;
    
    private double[] mu0_offset;
	
    //public ClusterAlgorithmOperator(MatrixParameter virusLocations, MatrixParameter mu, Parameter clusterLabels, Parameter K, double weight, Parameter virusOffsetsParameter, Parameter locationDrift_in, Parameter clusterOffsetsParameter) {
    public ClusterAlgorithmOperator(MatrixParameter virusLocations, MatrixParameter mu, Parameter clusterLabels, Parameter K, double weight, Parameter virusOffsetsParameter, Parameter clusterOffsetsParameter) {
    	
      	
    	System.out.println("Loading the constructor for ClusterAlgorithmOperator");
    	this.mu = mu;
    	this.K = K;
    	this.clusterLabels = clusterLabels;    	
    //	this.clusterLikelihood = clusterLikelihood;
        this.virusLocations = virusLocations;
        this.virusOffsetsParameter = virusOffsetsParameter;
    //    this.locationDrift = locationDrift_in;  //no longer need
        this.clusterOffsetsParameter = clusterOffsetsParameter;
        
        numdata = virusOffsetsParameter.getSize();
        System.out.println("numdata="+ numdata);
        
        
        int K_int = (int) K.getParameterValue(0);
        
        
        System.out.println("K_int=" + K_int);
        groupSize = new int[K_int];
        for(int i=0; i < K_int; i++){
        	groupSize[i] = 0;
        }
                
        
        for(int i=0; i < numdata; i++){
        	//System.out.println("i="+ i);
        	int index = (int) clusterLabels.getParameterValue(i);
        	groupSize[ index]++;
        }
    	
        for(int i=0; i < numdata;i++){
    		if(maxLabel < (int) clusterLabels.getParameterValue(i)){
    			maxLabel = (int) clusterLabels.getParameterValue(i);
    		}
    	}
        
        //NEED maxGROUP
        
        //for(int i=0; i < K_int; i++){
        	//System.out.println("groupSize=" + groupSize[i]);
        //}
        
        
        muLabels = new int[K_int];
        
        for(int i=0; i < maxLabel; i++){
        	int j=0;
            if(groupSize[i] >0){
            	muLabels[j] = i;
            	j++;
            }
        }
 
        //muLabels ...
        
        
        setWeight(weight);
        
        System.out.println("Finished loading the constructor for ClusterAlgorithmOperator");
   
    }
    

 
    /**
	 * change the parameter and return the log hastings ratio.
     */
    public final double doOperation() {
    	
    	//System.out.println("Do operation!");
    	
    	
  		double[] zeroVector2D = {0,0};
  		double[][] identityMatrix2D = new double[][]{
  			  { 1, 0 },
  			  { 0, 1 }};
  		
  		double[][] sigmaSqMatrix2D = new double[][]{
    			  { sigmaSq, 0 },
    			  { 0, sigmaSq }};
    	
    	double logHastingRatio = 0;  	
    	double chooseOperator = Math.random();
    	
    	int K_int = (int) K.getParameterValue(0);
    	
    	
    	
    	double[] original_groupSize = new double[groupSize.length];
		//recalculate groupSize
		for(int i=0; i < groupSize.length; i++){
			original_groupSize[i] = 0;
		}
		for(int i=0; i < numdata; i++){
    		int label =  (int) clusterLabels.getParameterValue(i);
    		original_groupSize[label  ]++;
		}

		
//		for(int i=0; i < K_int; i++){
//			System.out.println("group " + i + " has size=" + original_groupSize[i]);
//		}
		
  		
			for(int i=0; i < K_int; i++){
	      		double muk_0 = mu.getParameter(i).getParameterValue(0);
	      		double muk_1 = mu.getParameter(i).getParameterValue(1);  				
				//System.out.println("size=" +   groupSize[i] +   " mu_k_0=" + muk_0+ " , muk_1=" + muk_1);
			}
			
    	    //	System.out.println("propose a change in mu only");
			if(chooseOperator < 0.5){
		//	if(chooseOperator < 1){
	    		//change nothing
	    		isMoveMu = 1;
   	    		int groupSelect = (int) Math.floor( Math.random()* K_int );
   	    			groupSelectedChange = groupSelect;
  	    		int dimSelect = (int) Math.floor( Math.random()* 2 );   		
  	    			dimSelectChange = dimSelect;
	//    			 System.out.println("Group selected = " + groupSelectedChange + " mu=" + mu.getParameter(groupSelectedChange).getParameterValue(0)+ "\t" + mu.getParameter(groupSelectedChange).getParameterValue(1) + "    (before change...)"  );  	    		
	    		double change = Math.random()*2-1 ; 	
	    				//System.out.println(change);
      			double originalValue = mu.getParameter(groupSelect).getParameterValue(dimSelect);
      				originalValueChange = originalValue;
    			 mu.getParameter(groupSelect).setParameterValue(dimSelect, originalValue + change);
//   			 		 System.out.println("Group selected = " + groupSelectedChange + " mu=" + mu.getParameter(groupSelectedChange).getParameterValue(0)+ "\t" + mu.getParameter(groupSelectedChange).getParameterValue(1) + "    (propsed to...)");
	    		
	      		logHastingRatio = 0;
	    	}

			   // 		System.out.println("propose a change in both C and mu");
	    	else{		
	    		isMoveMu = 0;
	    			    		
	    		int virusIndex = (int) Math.floor( Math.random()*numdata );
	    		virusIndexChange = virusIndex;	    		
	    		int toBin = (int) Math.floor(Math.random()*K_int);
	    			//			System.out.println("toBin=" + toBin);
	    		int fromBin =  (int) clusterLabels.getParameterValue(virusIndex);
	    			// 		System.out.println("fromBin=" + fromBin);
	    	//	if(virusIndex < 5){
	    		//	System.out.println("virus " + virusIndex + "  from bin=" + fromBin + " to bin " + toBin);
	    	//	}
	    		clusterLabels.setParameterValue( virusIndex, toBin);   //the proposal

	    		
	    		
	    		
	    		
	    		
	    		//recalculate groupSize
	    		for(int i=0; i < groupSize.length; i++){
	    			groupSize[i] = 0;
	    		}
	    		for(int i=0; i < numdata; i++){
	        		int label =  (int) clusterLabels.getParameterValue(i);
	    			groupSize[label  ]++;
	    		}
	   		
	    		//special case that needs attention on the virus label
	    		if( (original_groupSize[fromBin] > 0) && ( groupSize[fromBin] == 0)){
	    			
	    			K.setParameterValue(0, K_int - 1);
	    			System.out.println("propose the fromBin " + fromBin + "becomes 0 in size - death of a bin");
	    			//actually that label is no longer used..
	    			double[] ranNormal =  MultivariateNormalDistribution.nextMultivariateNormalVariance( zeroVector2D, sigmaSqMatrix2D);
	    			
	    			mu.getParameter(fromBin).setParameterValue(0, ranNormal[0]);
	    			mu.getParameter(fromBin).setParameterValue(1, ranNormal[1]);
	    			//logHastingRatio += 0;  //this move doesn't change
	    		}
	    		
	    		//birth of a new bin.. assign an offset to it
	    		if( (original_groupSize[toBin] == 0) && (groupSize[toBin] == 1)){
	    			
	    			K.setParameterValue(0, K_int + 1);
	    			
	    			System.out.println("propose the birth of bin" + toBin);
	       			double offset = 0;
//	       			double drift = locationDrift.getParameterValue(0); // no longer need to do this here
	    //   			System.out.println("drift=" + drift);
	    	         if (virusOffsetsParameter != null) {
	    	            //	System.out.print("virus Offeset Parameter present"+ ": ");
	    	            //	System.out.print( virusOffsetsParameter.getParameterValue(i) + " ");
	    	            //	System.out.print(" drift= " + drift + " ");
	    	  //              offset = drift * virusOffsetsParameter.getParameterValue(virusIndex);
	    	                		//make sure that it is equivalent to double offset  = year[virusIndex] - firstYear;
	    	            }
	    	            else{
	    	            	System.out.println("virus Offeset Parameter NOT present. We expect one though. Something is wrong.");
	    	            }
	    			double[] ranNormal =  MultivariateNormalDistribution.nextMultivariateNormalVariance( zeroVector2D, sigmaSqMatrix2D);
	    			mu.getParameter(toBin).setParameterValue(0, ranNormal[0] ); // no need to assign offset anymore.. it's getting taken care of in the ClusterViruses by default
	    	//		mu.getParameter(toBin).setParameterValue(0, ranNormal[0] + offset);
	    			mu.getParameter(toBin).setParameterValue(1, ranNormal[1]);	    			
	    			
	    			//this move should change the Hasting Ratio!
	    			//CODE HERE
	    			
	    		}
	    		
	    	} //else    	 
			
	    	

	    /*   	
			for(int i=0; i < K_int; i++){
    			double muValue = mu.getParameter(i).getParameterValue(0);
    			double muValue2 = mu.getParameter(i).getParameterValue(1);
    			System.out.println("Group " + i + "\t" + muValue + "\t" + muValue2);
			}
			System.out.println("=============================");
			*/

			
			
			//change the mu in the toBin and fromBIn
			//borrow from getLogLikelihood:

			double[] meanYear = new double[K_int];
			double[] groupCount = new double[K_int];
			for(int i=0; i < numdata; i++){
				int label = (int) clusterLabels.getParameterValue(i);
				double year  = 0;
		        if (virusOffsetsParameter != null) {
		            //	System.out.print("virus Offeset Parameter present"+ ": ");
		            //	System.out.print( virusOffsetsParameter.getParameterValue(i) + " ");
		            //	System.out.print(" drift= " + drift + " ");
		                year = virusOffsetsParameter.getParameterValue(i);   //just want year[i]
		                		//make sure that it is equivalent to double offset  = year[virusIndex] - firstYear;
		            }
		            else{
		            	System.out.println("virus Offeset Parameter NOT present. We expect one though. Something is wrong.");
		            }
				meanYear[ label] = meanYear[ label] + year;
				
				groupCount[ label  ] = groupCount[ label ]  +1; 
			}
			int maxLabel=0;
			for(int i=0;i< numdata; i++){
				if(maxLabel < (int) clusterLabels.getParameterValue(i)){
					maxLabel = (int) clusterLabels.getParameterValue(i);
				}
			}
			
			for(int i=0; i <= maxLabel; i++){
				meanYear[i] = meanYear[i]/groupCount[i];
				//System.out.println(meanYear[i]);
			}
			
	
			//System.out.println("beta=" + beta);
			//beta = 1;

			mu0_offset = new double[maxLabel+1];
			//double[] mu1 = new double[maxLabel];
					
			
			//System.out.println("maxLabel=" + maxLabel);
			//now, change the mu..
			for(int i=0; i <= maxLabel; i++){
				//System.out.println(meanYear[i]*beta);
				mu0_offset[i] =  meanYear[i];
			//	System.out.println("group " + i + "\t" + mu0_offset[i]);
			}	
		//		System.out.println("=====================");
			
			
						

			//Set  the vLoc to be the corresponding mu values , and clusterOffsetsParameter to be the corresponding offsets
	    	//virus in the same cluster has the same position
	    	for(int i=0; i < numdata; i++){
	        	int label = (int) clusterLabels.getParameterValue(i);
	    		Parameter vLoc = virusLocations.getParameter(i);
	    		
	    		//setting the virus locs to be equal to the corresponding mu
	    			double muValue = mu.getParameter(label).getParameterValue(0);
	    			vLoc.setParameterValue(0, muValue);

	    		  double	muValue2 = mu.getParameter(label).getParameterValue(1);
	   				vLoc.setParameterValue(1, muValue2);
		
	   			
	   			//if we want to apply the mean year virus cluster offset to the cluster
	   			if(clusterOffsetsParameter != null){
	   			//setting the clusterOffsets to be equal to the mean year of the virus cluster
	   				// by doing this, the virus changes cluster AND updates the offset simultaneously
	   				clusterOffsetsParameter.setParameterValue( i , mu0_offset[label]);
	   			}
     					//	System.out.println("mu0_offset[label]=" + mu0_offset[label]);
     				//	System.out.println("clusterOffsets now becomes =" + clusterOffsetsParameter.getParameterValue(i) );   			
	    	}
	   // 	System.out.println("");
	    	
	    	
	    	//Hasting's Ratio is p(old |new)/ p(new|old)

	    	//System.out.println("Done doing operation!");
			
	    		    	
	    	
    	//return(logHastingRatio); //log hasting ratio
    	return(logHastingRatio);
    	
    }
    	
    	
    public void accept(double deviation) {
    	super.accept(deviation);

    	/*
    	if(isMoveMu==1){
    		numAcceptMoveMu++;
    		numProposeMoveMu++;
        	System.out.println("% accept move Mu = " + numAcceptMoveMu/(double)numProposeMoveMu);
    	}
    	else{    	   
    		numAcceptMoveC++;
    		numProposeMoveC++;
        	System.out.println("% accept move C = " + numAcceptMoveC/(double)numProposeMoveC);
    	}
    	*/  
    	        	
    	//	if(virusIndexChange <5){
    //		System.out.println("     -  Accepted!");
    	//	}
          	
    }
    
    public void reject(){
    	super.reject();
 
    	
    	/*
    	//manually change mu back..
    	if(isMoveMu==1){
			 mu.getParameter(groupSelectedChange).setParameterValue(dimSelectChange, originalValueChange);
    	}
    	//manually change all the affected vLoc back...
    	for(int i=0; i < numdata; i++){
        	int label = (int) clusterLabels.getParameterValue(i);
    		Parameter vLoc = virusLocations.getParameter(i);   		
    		//	double muValue = mu.getParameter(label).getParameterValue(0);
    		//	vLoc.setParameterValue(0, muValue);
    		//  double	muValue2 = mu.getParameter(label).getParameterValue(1);
   			//	vLoc.setParameterValue(1, muValue2);
	
 			clusterOffsetsParameter.setParameterValue( i , mu0_offset[label]);   			
    	}
    	*/
    	
    	
    	/*
    	if(isMoveMu==1){
    		numProposeMoveMu++;
        	System.out.println("% accept move Mu = " + numAcceptMoveMu/(double)numProposeMoveMu);
    	}
    	else{    	   
    		numProposeMoveC++;
        	System.out.println("% accept move C = " + numAcceptMoveC/(double)numProposeMoveC);
    	}
    	*/
    	//if(virusIndexChange < 5){
//		System.out.println("     -      Rejected!");
    	//}
      	
      	
      	/*
      	for(int i=0; i < numdata; i++){
      		Parameter vLoc = virusLocations.getParameter(i);

      		if( vLoc.getParameterValue(0) != old_vLoc0[i]){

      			System.out.println("virus " + i + " is different: " + vLoc.getParameterValue(0) + " and " + old_vLoc0[i]);
      		}
      		
      		//System.out.println(old_vLoc0[i] + ", " + old_vLoc1[i]);
      		vLoc.setParameterValue(0, old_vLoc0[i]);
      		vLoc.setParameterValue(1, old_vLoc1[i]);
      		
		}
      	*/
  		//System.exit(0);

      	
      	

    }
    
	

           public final static String CLUSTERALGORITHM_OPERATOR = "ClusterAlgorithmOperator";

              
            //MCMCOperator INTERFACE
            public final String getOperatorName() {
                return CLUSTERALGORITHM_OPERATOR;
            }

            public final void optimize(double targetProb) {

                throw new RuntimeException("This operator cannot be optimized!");
            }

            public boolean isOptimizing() {
                return false;
            }

            public void setOptimizing(boolean opt) {
                throw new RuntimeException("This operator cannot be optimized!");
            }

            public double getMinimumAcceptanceLevel() {
                return 0.1;
            }

            public double getMaximumAcceptanceLevel() {
                return 0.4;
            }

            public double getMinimumGoodAcceptanceLevel() {
                return 0.20;
            }

            public double getMaximumGoodAcceptanceLevel() {
                return 0.30;
            }

            public String getPerformanceSuggestion() {
                if (Utils.getAcceptanceProbability(this) < getMinimumAcceptanceLevel()) {
                    return "";
                } else if (Utils.getAcceptanceProbability(this) > getMaximumAcceptanceLevel()) {
                    return "";
                } else {
                    return "";
                }
            }

        
           
        

            public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {
            	

                public final static String VIRUSLOCATIONS = "virusLocations";
            	public final static String  MU = "mu";
            	public final static String CLUSTERLABELS = "clusterLabels";
            	public final static String K = "k";
            	public final static String OFFSETS = "offsets";
      //     	public final static String LOCATION_DRIFT = "locationDrift"; //no longer need
            	
                public final static String CLUSTER_OFFSETS = "clusterOffsetsParameter";


                public String getParserName() {
                    return CLUSTERALGORITHM_OPERATOR;
                }

                /* (non-Javadoc)
                 * @see dr.xml.AbstractXMLObjectParser#parseXMLObject(dr.xml.XMLObject)
                 */
                public Object parseXMLObject(XMLObject xo) throws XMLParseException {
                	
                	//System.out.println("Parser run. Exit now");
                	//System.exit(0);

                    double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);

                    
                    XMLObject cxo = xo.getChild(VIRUSLOCATIONS);
                        MatrixParameter virusLocations = (MatrixParameter) cxo.getChild(MatrixParameter.class);

                        
                        cxo = xo.getChild(MU);
                        MatrixParameter mu = (MatrixParameter) cxo.getChild(MatrixParameter.class);

                        cxo = xo.getChild(CLUSTERLABELS);
                        Parameter clusterLabels = (Parameter) cxo.getChild(Parameter.class);

                        cxo = xo.getChild(K);
                        Parameter k = (Parameter) cxo.getChild(Parameter.class);
                        
                        cxo = xo.getChild(OFFSETS);
                        Parameter offsets = (Parameter) cxo.getChild(Parameter.class);
 
//                        cxo = xo.getChild(LOCATION_DRIFT);
//                        Parameter locationDrift = (Parameter) cxo.getChild(Parameter.class);
                        
                        Parameter clusterOffsetsParameter = null;
                        if (xo.hasChildNamed(CLUSTER_OFFSETS)) {
                        	clusterOffsetsParameter = (Parameter) xo.getElementFirstChild(CLUSTER_OFFSETS);
                        }

                        
                        
                    //return new ClusterAlgorithmOperator(virusLocations, mu, clusterLabels, k, weight, offsets, locationDrift, clusterOffsetsParameter);
                        return new ClusterAlgorithmOperator(virusLocations, mu, clusterLabels, k, weight, offsets,  clusterOffsetsParameter);

                }

                //************************************************************************
                // AbstractXMLObjectParser implementation
                //************************************************************************

                public String getParserDescription() {
                    return "An operator that picks a new allocation of an item to a cluster under the Dirichlet process.";
                }

                public Class getReturnType() {
                    return ClusterAlgorithmOperator.class;
                }


                public XMLSyntaxRule[] getSyntaxRules() {
                    return rules;
                }

                private final XMLSyntaxRule[] rules = {
                        AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
                        
                        new ElementRule(VIRUSLOCATIONS, Parameter.class),
                        new ElementRule(MU, Parameter.class),
                        new ElementRule(CLUSTERLABELS, Parameter.class),
                        new ElementRule(K, Parameter.class),
                        new ElementRule(OFFSETS, Parameter.class),
                  //      new ElementRule(LOCATION_DRIFT, Parameter.class), //no longer needed
   //                    
                       new ElementRule(CLUSTER_OFFSETS, Parameter.class, "Parameter of cluster offsets of all virus"),  // no longer REQUIRED
                
        
            };
            
            };


        
            public int getStepCount() {
                return 1;
            }

        }





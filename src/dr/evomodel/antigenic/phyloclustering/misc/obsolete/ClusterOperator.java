package dr.evomodel.antigenic.phyloclustering.misc.obsolete;

import dr.evomodel.antigenic.NPAntigenicLikelihood;
import dr.inference.model.Likelihood;
import dr.inference.model.MatrixParameter;
import dr.inference.model.Parameter;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;
import dr.xml.*;



public class ClusterOperator extends SimpleMCMCOperator  {



            public final static String CLUSTER_OPERATOR = "clusterOperator";
            public static final String WINDOW_SIZE = "windowSize";
            
    
            public NPAntigenicLikelihood modelLikelihood;           
            private Parameter assignments = null;
            private Parameter links = null;
            
            private MatrixParameter virusLocations = null;
            private double windowSize = 1;

                    
            
            public ClusterOperator(Parameter assignments, MatrixParameter virusLocations, double weight, double windowSize, NPAntigenicLikelihood Likelihood, Parameter links) {
            	//System.out.println("Constructor");
                this.links = links;           	
                this.assignments = assignments;
                this.virusLocations = virusLocations;
                this.windowSize = windowSize;
                this.modelLikelihood = Likelihood;
                
                setWeight(weight);
                
                
                //load in the virusLocations
                System.out.println((int)assignments.getParameterValue(0));
                Parameter vLoc = virusLocations.getParameter(0);
        
        		//double d0=  vLoc.getParameterValue(0);
        		//double d1= vLoc.getParameterValue(1);
        		//System.out.println( d0 + "," + d1);
        		
                
                //System.exit(0);

			}




            //assume no boundary
            /**
             * change the parameter and return the hastings ratio.
             */
            public final double doOperation() {

                double moveWholeClusterOrChangeOnePoint = MathUtils.nextDouble();
            	            	
            	//double moveWholeClusterProb = 0.2;
            	double moveWholeClusterProb = 0.5;
            	
            	int numSamples = assignments.getDimension();           			

//           	if(moveWholeClusterOrChangeOnePoint <moveWholeClusterProb){
      //      	System.out.print("Group changes: ");	
            	//System.out.println("Move a group of points in a cluster together");
            	int target = MathUtils.nextInt(numSamples);
            	int targetCluster = (int) assignments.getParameterValue(target);

 
            	//for(int i=0; i < numSamples; i++){
                   //Parameter vLoc = virusLocations.getParameter(i);
            		//double d0=  vLoc.getParameterValue(0);
            		//double d1= vLoc.getParameterValue(1);
            		//System.out.println( d0 + "," + d1);
            	//}
 
                // a random dimension to perturb
                int ranDim= MathUtils.nextInt(virusLocations.getParameter(0).getDimension());
                // a random point around old value within windowSize * 2
                double draw = (2.0 * MathUtils.nextDouble() - 1.0) * windowSize;
                
                
                
                for(int i=0; i < numSamples; i++){
                	if((int) assignments.getParameterValue(i)== targetCluster ){
                		//System.out.print("update "+ i +" from ");
                		Parameter vLoc = virusLocations.getParameter(i);
                		//System.out.print(vLoc.getParameterValue(index));
                		double newValue = vLoc.getParameterValue(ranDim) + draw;
                		vLoc.setParameterValue(ranDim, newValue);
                		//System.out.println(" to " + vLoc.getParameterValue(index));
                	}
                }
               
  //          }
            	
    /*        	
            	//change cluster assignment!
            	else{
            		System.out.print("Cluster assignment & Group move:");
            		//System.out.println("Change a single point to another cluster");
            		//moving the point without coupling the point is useless!                 
                    
            		
            		
            		//randomly pick a point to change cluster assignment
            		//move the point and couple it with change in cluster:
            		
                    int index = MathUtils.nextInt(numSamples);   		
                    int oldGroup = (int)assignments.getParameterValue(index);
                   
                  // Set index customer link to index and all connected to it to a new assignment (min value empty)                  
                    int minEmp = minEmpty(modelLikelihood.getLogLikelihoodsVector());
                    links.setParameterValue(index, index);
                    int[] visited = connected(index, links);
                    
                    int ii = 0;
                    while (visited[ii]!=0){
                            assignments.setParameterValue(visited[ii]-1, minEmp);
                            ii++;
                    }

                   //Adjust likvector for group separated
                   modelLikelihood.setLogLikelihoodsVector(oldGroup,modelLikelihood.getLogLikGroup(oldGroup) );
                   modelLikelihood.setLogLikelihoodsVector(minEmp,modelLikelihood.getLogLikGroup(minEmp) );
                   int maxFull = maxFull( modelLikelihood.getLogLikelihoodsVector());

                   
                   //pick new link
                   int k = MathUtils.nextInt(numSamples); //wonder if there is a way to make this more efficient at choosing good moves             
                   
                   links.setParameterValue(index, k);                  
                   int newGroup = (int)assignments.getParameterValue(k);
                   int countNumChanged = 0;
                   ii = 0;
                   while (visited[ii]!=0){
                           assignments.setParameterValue(visited[ii]-1, newGroup);
                           ii++;
                           countNumChanged++;
                   }
                   
                   System.out.print("changed="+countNumChanged+ " ");
                   
                   int targetCluster = (int) assignments.getParameterValue(k);                  
                   
                  
                   // updating conditional likelihood vector 
                   modelLikelihood.setLogLikelihoodsVector(newGroup, modelLikelihood.getLogLikGroup(newGroup));
                   if (newGroup!=minEmp){
                            modelLikelihood.setLogLikelihoodsVector(minEmp, 0);
                          
                   }

                         
            		
            		
            	//change the location
               // a random dimension to perturb
               int ranDim= MathUtils.nextInt(virusLocations.getParameter(0).getDimension());
               
               // a random point around old value within windowSize * 2
               double draw = (2.0 * MathUtils.nextDouble() - 1.0) * windowSize;
               
               
               int countNumLocChanged = 0;
               for(int i=0; i < numSamples; i++){
               	if((int) assignments.getParameterValue(i)== targetCluster ){
               		//System.out.print("update "+ i +" from ");
               		Parameter vLoc = virusLocations.getParameter(i);
               		//System.out.print(vLoc.getParameterValue(index));
               		double newValue = vLoc.getParameterValue(ranDim) + draw;
               		vLoc.setParameterValue(ranDim, newValue);
               		//System.out.println(" to " + vLoc.getParameterValue(index));
               		countNumLocChanged++;
               	}
               }
               System.out.print("loc changed="+countNumLocChanged+ " ");
 
            		
              // sampleMeans(maxFull); // I decide not to resample the means here.
           		
            }
            //Want to check that virusLocations WON't get updated if rejected
                
                
//            	System.out.println(target);
 //           	System.out.println(targetCluster);
            	//for(int i=0; i < assignments.getDimension();)
      */      	
            	
                return 0.0;
        
            }


            
            /*
             * find min Empty
             */
            
            public int minEmpty(double[] logLikVector){
            int isEmpty=0;
            int i =0;
                while (isEmpty==0){
            if(logLikVector[i]==0){
                isEmpty=1;}
            else { 
                if(i==logLikVector.length-1){isEmpty=1;}
                   i++;}
                }
            return i;
            }
           
        
            
            /*
             * find max Full
             */
            
            
            public int maxFull(double[] logLikVector){
                    int isEmpty=1;
                    int i =logLikVector.length-1;
                        while (isEmpty==1){
                    if(logLikVector[i]!=0){
                        isEmpty=0;}
                    else {i--;}
                        }
                    return i;
                    }
         /*
          * find customers connected to i
          */
            
            public int[] connected(int i, Parameter clusteringParameter){
                int n =  clusteringParameter.getDimension();
                int[] visited = new int[n+1]; 
                    visited[0]=i+1;
                    int tv=1;
                    
                    for(int j=0;j<n;j++){
                        if(visited[j]!=0){
                                int curr = visited[j]-1;
                                
                                /*look forward
                                */
                                
                        int forward = (int) clusteringParameter.getParameterValue(curr);
                        visited[tv] = forward+1;
                        tv++;
                            // Check to see if is isn't already on the list

                        for(int ii=0; ii<tv-1; ii++){
                        if(visited[ii]==forward+1){
                                tv--;
                                visited[tv]=0;
                        }
                        }


                        /*look back
                        */
                        for (int jj=0; jj<n;jj++){
                                if((int)clusteringParameter.getParameterValue(jj)==curr){ 
                                        visited[tv]= jj+1;
                                        tv++;
                                        
                                        for(int ii=0; ii<tv-1; ii++){
                                        if(visited[ii]==jj+1){
                                                tv--;
                                                visited[tv]=0;
                                        }               
                                }
          
                        }
                        }
                        
                        }}
                return visited;
                
            }
                        
            
            
            public void accept(double deviation) {
            	super.accept(deviation);
  //          	System.out.println("Accepted!");
            	
            }
            
            public void reject(){
            	super.reject();
    //        	System.out.println("Rejected");
            }
            
            
            
            //MCMCOperator INTERFACE
            public final String getOperatorName() {
                return CLUSTER_OPERATOR;
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
                public final static String ASSIGNMENTS = "assignments";
                public final static String VIRUS_LOCATIONS = "virusLocations";

                public final static String LIKELIHOOD = "likelihood";
                public final static String LINKS = "links";
                
                
                public String getParserName() {
                    return CLUSTER_OPERATOR;
                }

                /* (non-Javadoc)
                 * @see dr.xml.AbstractXMLObjectParser#parseXMLObject(dr.xml.XMLObject)
                 */
                public Object parseXMLObject(XMLObject xo) throws XMLParseException {

                    double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
              //      double windowSize = xo.getDoubleAttribute(WINDOW_SIZE);

                    
                    XMLObject cxo =  xo.getChild(ASSIGNMENTS);
                        Parameter assignments = (Parameter) cxo.getChild(Parameter.class);  
                        
                        MatrixParameter virusLocationsParameter = null;
                        if (xo.hasChildNamed(VIRUS_LOCATIONS)) {
                            virusLocationsParameter = (MatrixParameter) xo.getElementFirstChild(VIRUS_LOCATIONS);
                        }
                   
                        cxo = xo.getChild(LINKS);
                        Parameter links = (Parameter) cxo.getChild(Parameter.class);

                        cxo = xo.getChild(LIKELIHOOD);
                        NPAntigenicLikelihood likelihood = (NPAntigenicLikelihood)cxo.getChild(NPAntigenicLikelihood.class);
                                                
                        //set window size to be 1
                    return new ClusterOperator(assignments, virusLocationsParameter, weight, 1.0, likelihood, links);
          //          public ClusterOperator(Parameter assignments, MatrixParameter virusLocations, double weight, double windowSize, NPAntigenicLikelihood Likelihood, Parameter links) {


                }

                //************************************************************************
                // AbstractXMLObjectParser implementation
                //************************************************************************

                public String getParserDescription() {
                    return "An operator that picks a new allocation of an item to a cluster under the Dirichlet process.";
                }

                public Class getReturnType() {
                    return ClusterOperator.class;
                }


                public XMLSyntaxRule[] getSyntaxRules() {
                    return rules;
                }

                private final XMLSyntaxRule[] rules = {
                        AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
                   //     AttributeRule.newDoubleRule(WINDOW_SIZE),                        
                        new ElementRule(ASSIGNMENTS,
                                new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
                       new ElementRule(VIRUS_LOCATIONS, MatrixParameter.class, "Parameter of locations of all virus"),
                       new ElementRule(LINKS,
                               new XMLSyntaxRule[]{new ElementRule(Parameter.class)}), 
                        new ElementRule(LIKELIHOOD, new XMLSyntaxRule[] {
                                       new ElementRule(Likelihood.class),}, true),                              
        
            };
            
            };            
            


       
            public int getStepCount() {
                return 1;
            }

        }





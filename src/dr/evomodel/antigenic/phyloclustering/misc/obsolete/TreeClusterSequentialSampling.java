
package dr.evomodel.antigenic.phyloclustering.misc.obsolete;

import java.util.HashMap;
import java.util.LinkedList;
import java.util.Map;

import dr.evolution.tree.NodeRef;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.MatrixParameter;
import dr.inference.model.Parameter;
import dr.inference.operators.GibbsOperator;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;
import dr.xml.*;


/**
 * A Gibbs operator for allocation of items to clusters under a distance dependent Chinese restaurant process.
 *
 * @author Charles Cheung
 * @author Trevor Bedford
 */

public class TreeClusterSequentialSampling extends SimpleMCMCOperator implements GibbsOperator{
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
    private AGLikelihoodTreeCluster clusterLikelihood = null;
   
    private int groupSelectedChange = -1;
    private int virusIndexChange = -1;
    private double originalValueChange = -1;
    private int dimSelectChange = -1;
    
    private double[] mu0_offset;
    
    private Parameter breakPoints = null;

    private int binSize=20;
 
    
    private Parameter status;
    
    private TreeModel treeModel;

    int []membershipToClusterLabelIndexes = null;
    
    private int numNodes;
    
    
   // private int[] piIndicator = new int[numSites];

	
    //public ClusterAlgorithmOperator(MatrixParameter virusLocations, MatrixParameter mu, Parameter clusterLabels, Parameter K, double weight, Parameter virusOffsetsParameter, Parameter locationDrift_in, Parameter clusterOffsetsParameter) {
    public TreeClusterSequentialSampling(MatrixParameter virusLocations, MatrixParameter mu, Parameter clusterLabels, Parameter K, double weight, Parameter virusOffsetsParameter, Parameter clusterOffsetsParameter, Parameter breakPointsParameter, Parameter statusParameter, TreeModel treeModel_in, AGLikelihoodTreeCluster clusterLikelihood_in) {
    	
      	
    	System.out.println("Loading the constructor for Sequential sampler");
    	
    	this.clusterLikelihood = clusterLikelihood_in;
		this.treeModel= treeModel_in;
    	this.mu = mu;
    	this.K = K;
    	this.clusterLabels = clusterLabels;    	
    //	this.clusterLikelihood = clusterLikelihood;
        this.virusLocations = virusLocations;
        this.virusOffsetsParameter = virusOffsetsParameter;
    //    this.locationDrift = locationDrift_in;  //no longer need
        this.clusterOffsetsParameter = clusterOffsetsParameter;
    	this.breakPoints= breakPointsParameter;   
    	
    	this.status = statusParameter;

        
        numdata = virusOffsetsParameter.getSize();
        System.out.println("numdata="+ numdata);
        
        
        int K_int = (int) K.getParameterValue(0);
        
        
        System.out.println("K_int=" + K_int);
        groupSize = new int[binSize];
        for(int i=0; i < binSize; i++){
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
        
        
        muLabels = new int[binSize];
        
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
    	
    	System.out.println("do operation of sequential sampling");
    	
    	setMembershipToClusterLabelIndexes(); //run once in case the tree changes.
    	
    	numNodes = treeModel.getNodeCount(); 
    	
    	updateK();
    	int K_int = (int) K.getParameterValue(0); 
    	
    	//System.out.println("K_int is " + K_int);
    	
    	int[] cutNode = new int[K_int]; //to store which nodes to add...
   
        int[] oldclusterLabelArray = new int[numNodes]; //everything 0.
        int[] clusterLabelArray = null;
        
        int[] onPoints = new int[K_int];
        int numOn=0;
        for(int i=0; i < binSize; i++){
        	if( (int)status.getParameterValue(i) ==1 ){
        		onPoints[numOn] = i;
        		numOn++;
        	}
        }
        
        
        //Verify relationship 
        //P(H|E0=a, Y, mu) = sum P(H|E0=a, E1=e1, Y, mu) x P(E1=e1| E0=a)
        

        cutNode = new int[1];
    	cutNode[0] = 785;
    	onPoints = new int[1];
    	onPoints[0] = 0;
       	//now translate the cutNodes into the breakpoints that are on [although here we don't use the indicators - elsewhere does
    	resetStatusAndBreakpointsGivenCutNodes(cutNode, onPoints);//need to reset the status breakpoints by the testCutNode
 
    	
    	//use the tree to re-partition according to the change.
		clusterLabelArray = setClusterLabelsByTestCutNodeByNodeOrder(cutNode);
		relabelClusterLabels(clusterLabelArray, oldclusterLabelArray); //will move it out
		for(int i=0; i < numdata; i++){
			clusterLabels.setParameterValue(i, clusterLabelArray[i]); 	//set cluster label parameter 					
		}
		
		//setVirusLocationAndOffsets(); //set virus locations, given the clusterlabels parameter
		setVirusLocationAutoCorrelatedModel();		
		
		double topLogProb = clusterLikelihood.getLogLikelihood();
		System.out.println(topLogProb);
		   		
    			

		//now test the sum...
        cutNode = new int[2];
    	cutNode[0] = 785;
    	onPoints = new int[2];
    	onPoints[0] = 0;
    	onPoints[1] = 1;
 
		double []logNumerator = new double[numNodes];
		//calculate the distribution for calculating introducing an excision point in each node
		for(int curTest=0; curTest < numNodes; curTest++){	
			
			int hasBeenAdded = checkSiteHasBeenAddedToOnIndicators(curTest); 			//check if a site has already been added
  			if(hasBeenAdded ==0){
  				
  				cutNode = new int[2];
  				cutNode[0] = 785;
  				cutNode[1] = curTest;

  		    	resetStatusAndBreakpointsGivenCutNodes(cutNode, onPoints);//need to reset the status breakpoints by the testCutNode
  		    	//use the tree to re-partition according to the change.
  				clusterLabelArray = setClusterLabelsByTestCutNodeByNodeOrder(cutNode); //note that instead of using the indicators, it uses the testCutNodes directly
  				relabelClusterLabels(clusterLabelArray, oldclusterLabelArray); //will move it out
  				
					//set cluster label parameter for testing 					
  				for(int i=0; i < numdata; i++){
  					clusterLabels.setParameterValue(i, clusterLabelArray[i]);
  				}
		  				  				
   				//setVirusLocationAndOffsets();  //this uses the clusterLabels parameter
   				setVirusLocationAutoCorrelatedModel(); //which depends on the status and breakpoints
   				
   				   							    		  					  				
   				logNumerator[curTest] += clusterLikelihood.getLogLikelihood() ; 	//Calculate likelihood
  			}
  			else{
	  			logNumerator[curTest] = Double.NEGATIVE_INFINITY;
  				System.out.println("Don't calculate for node" + cutNode[0]);
  			}

		} //finished curTest
		
		
		double answer = verifyAssumption(topLogProb, logNumerator);
		
		
		System.out.println("the ratio is " + answer);
		

  		    	
  		    	
  		    	
  		    	
		
		
    	
    	System.exit(0);

        
        
        
        
        
        
        
        
        //switch
        
        double tmpMu1 = mu.getParameter(onPoints[7] +1).getParameterValue(0);
        double tmpMu2 = mu.getParameter(onPoints[7] +1).getParameterValue(1);
        
        mu.getParameter(onPoints[7] +1).setParameterValue(0, mu.getParameter(onPoints[8] +1).getParameterValue(0));
        mu.getParameter(onPoints[7] +1).setParameterValue(1, mu.getParameter(onPoints[8] +1).getParameterValue(1));
        
        mu.getParameter(onPoints[8] +1).setParameterValue(0, tmpMu1);
        mu.getParameter(onPoints[8] +1).setParameterValue(1, tmpMu2);
        
        
    	for(int curNode=8; curNode < K_int; curNode++){
    		
    		cutNode[0] = 785;
    		cutNode[1] = 775;
    		cutNode[2] = 763;
    		cutNode[3] = 697;
    		cutNode[4] = 747;
    		cutNode[5] = 679;
    		cutNode[6] = 662;
    		//cutNode[7] = 521;  //although 526 will be better than 521
    		cutNode[7] = 638;

    		
    		
    		//calculate the conditional distribution of the curNode, given the current set of nodes
    		double []logNumeratorProb = new double[numNodes];
    		
    		//calculate the distribution for calculating introducing an excision point in each node
    		for(int curTest=0; curTest < numNodes; curTest++){	
    			
    			int hasBeenAdded = checkSiteHasBeenAddedToOnIndicators(curTest); 			//check if a site has already been added
	  			if(hasBeenAdded ==0){
	  	  		  	int[] testCutNode = new int[curNode+1];	    			//create the testCutNode, with adding the current test node
	  	  		  	for(int i=0; i < curNode; i++){
	  	  		  		testCutNode[i] = cutNode[i];
	  	  		  	}
	  	    		testCutNode[curNode] = curTest;

	  		    	resetStatusAndBreakpointsGivenCutNodes(testCutNode, onPoints);//need to reset the status breakpoints by the testCutNode

	  		    	
	  		    	/*
	  		    	System.out.print("Currently selected: [");
	  		    	for(int i=0; i < binSize; i++){
	  		    		if((int) status.getParameterValue(i) ==1){
	  		    			System.out.print(  (int) breakPoints.getParameterValue(i) + ",");
	  		    		}
	  		    	}
	  		    	System.out.println("]");
	  		    	*/
	  	    		//I suspect I need to change this as I modify the code..
	  		       	//set the indicators, based on the cutnodes
	  		    	
	  		    	/*
	  		    	int addCount=0;
	  		    	for(int i=0; i < binSize; i++){
	  		    		if( (int)status.getParameterValue(i) == 1){
	  		    			breakPoints.setParameterValue(i, testCutNode[addCount]);
	  		    			addCount++;
	  		    		}
	  		    		if(addCount == (curNode+1)){
	  		    			break;
	  		    		}
	  		    	}
	  		    	*/
	  		    	//System.out.println("currently added " + addCount + " nodes");

	  		    	
	  	    		
	  		    	//use the tree to re-partition according to the change.
	  				clusterLabelArray = setClusterLabelsByTestCutNodeByNodeOrder(testCutNode); //note that instead of using the indicators, it uses the testCutNodes directly
	  				relabelClusterLabels(clusterLabelArray, oldclusterLabelArray); //will move it out
	  				
  					//set cluster label parameter for testing 					
	  				for(int i=0; i < numdata; i++){
	  					clusterLabels.setParameterValue(i, clusterLabelArray[i]);
	  				}
	  					  				
	  				
	   				//setVirusLocationAndOffsets();  //this uses the clusterLabels parameter
	   				setVirusLocationAutoCorrelatedModel(); //which depends on the status and breakpoints
	   				
	   				
	   				if(curNode == 0 && curTest == 0){
			    		for(int i=0; i < numdata; i++){
			    			Parameter v = virusLocations.getParameter(i);
			    			//v.setParameterValue(0, 0);
			    			//v.setParameterValue(1, 0);
			    		}
	   				}
	   							    		  					  				
   			    	logNumeratorProb[curTest] = clusterLikelihood.getLogLikelihood(); 	//Calculate likelihood
   			    	
   			    	//if(curNode == 0 && curTest == 0){
   			    		
   			    		
   			    		//for(int i=0; i < numdata; i++){
   			    		//	Parameter v = virusLocations.getParameter(i);
   			    		//	System.out.print(v.getParameterValue(0) +"," + v.getParameterValue(1)+"\t");
   			    		//}
   			    		//System.out.println("");
   			    		
   			    		//System.out.println(" * " + logNumeratorProb[curTest]);
   			    		//System.exit(0);
   			    	//}
	  			}
	  			else{
		  			logNumeratorProb[curTest]  = Double.NEGATIVE_INFINITY; //dummy probability
	  			}

    		} //finished curTest
    		
    		 double []condDistribution = calculateConditionalProbabilityGivenLogNumeratorProb(logNumeratorProb);
    		
    		 
    		 for(int i=0; i < numNodes; i++){
    			// if(condDistribution[i] > 0.0000001){
    				 System.out.println("node " + i + " p=" + condDistribution[i]);
    			// }
    		 }
    		 System.out.println("===============================");
    		 
    		 
    		// System.exit(0);
    		 
    		 int site_add = MathUtils.randomChoicePDF(condDistribution); //sample a site, given the conditioanl distribution
    		
  		  	//update the cutNode using a temporary array newCutNode - first copy the existing element, then add the new site
  		  	int[] newCutNode = new int[curNode+1];
  		  	for(int i=0; i < curNode; i++){
  		  		newCutNode[i] = cutNode[i];
  		  	}
    		newCutNode[curNode] = site_add; //add the new site.
    		cutNode = newCutNode;
    								
						
			//now , after adding the node to cutNode, we need to update clusterLabel such that the labeling is consistent with the old cluster labels.. 
				// because ie. cluster i always get mu i.
			clusterLabelArray = setClusterLabelsByTestCutNodeByNodeOrder(newCutNode);
			relabelClusterLabels(clusterLabelArray, oldclusterLabelArray); //will move it out		
			
			oldclusterLabelArray = clusterLabelArray; // keep the oldcluster label to build on top of it..			
    	} //curNode    	

    	
    	
    	
    	//====================================================================================================
		//After finishing the proposal
		//====================================================================================================

    	//Display:
    	printCutNode(cutNode);

       	//now translate the cutNodes into the breakpoints that are on [although here we don't use the indicators - elsewhere does
    	resetStatusAndBreakpointsGivenCutNodes(cutNode, onPoints);//need to reset the status breakpoints by the testCutNode
    	/*
    	int addCount=0;
    	for(int i=0; i < binSize; i++){
    		if( (int)status.getParameterValue(i) == 1){
    			breakPoints.setParameterValue(i, cutNode[addCount]);
    			addCount++;
    		}
    	}
    	*/
 
    	
    	//use the tree to re-partition according to the change.
		clusterLabelArray = setClusterLabelsByTestCutNodeByNodeOrder(cutNode);
		relabelClusterLabels(clusterLabelArray, oldclusterLabelArray); //will move it out
		for(int i=0; i < numdata; i++){
			clusterLabels.setParameterValue(i, clusterLabelArray[i]); 	//set cluster label parameter 					
		}
		
		//setVirusLocationAndOffsets(); //set virus locations, given the clusterlabels parameter
		setVirusLocationAutoCorrelatedModel();
		
		
		
		System.out.println(clusterLikelihood.getLogLikelihood());
		
		
		
		
		
		
		
		
		//add manually to test...
		cutNode[0] = 785;
		cutNode[1] = 775;
		cutNode[2] = 763;
		cutNode[3] = 697;
		cutNode[4] = 747;
		cutNode[5] = 679;
		cutNode[6] = 662;
		//cutNode[7] = 521;
		//cutNode[8] = 638;
		
		cutNode[7] = 638;
		cutNode[8] = 521;
		
    	resetStatusAndBreakpointsGivenCutNodes(cutNode, onPoints);//need to reset the status breakpoints by the testCutNode
    	//use the tree to re-partition according to the change.
		clusterLabelArray = setClusterLabelsByTestCutNodeByNodeOrder(cutNode);
		relabelClusterLabels(clusterLabelArray, oldclusterLabelArray); //will move it out
		for(int i=0; i < numdata; i++){
			clusterLabels.setParameterValue(i, clusterLabelArray[i]); 	//set cluster label parameter 					
		}
		
		//setVirusLocationAndOffsets(); //set virus locations, given the clusterlabels parameter
		setVirusLocationAutoCorrelatedModel();
	
		System.out.println(clusterLikelihood.getLogLikelihood());
		System.exit(0);
		
		
		
		
		//NEED TO SET THE OTHER BREAKPOINTS TO MAKE SURE THERE ARE NO DUPLICATES AT THE END...
		//SO OTHER PROPOSALS WOULD FUNCTION OK OK
		
 
		
		return(Double.POSITIVE_INFINITY); //it should be anything... always accept for the Gibbs move.    	
    }
    	
   





	private void resetStatusAndBreakpointsGivenCutNodes(int[] testCutNode, int[] onPoints) {
		
		for(int i=0; i < binSize; i++){
			status.setParameterValue(i, 0);
			breakPoints.setParameterValue(i, -1);
		}
		
    	int numOn = testCutNode.length;
    	int countOn=0;
    	for(int i=0; i < binSize; i++){
   			if(countOn < numOn){
   				status.setParameterValue(onPoints[countOn], 1);
    			breakPoints.setParameterValue(onPoints[countOn],testCutNode[countOn]); //reset breakPoints accordingly
    			countOn++;
    		}
    		
    	}
	}



	private void updateK() {

    	//K is changed accordingly..
		int K_count = 0; //K_int gets updated
		for(int i=0; i < binSize; i++){
			K_count += (int) status.getParameterValue(i);
		}
		//System.out.println("K now becomes " + K_count);
		K.setParameterValue(0, K_count); //update 
 			
				
	}

	
	private void printCutNode(int[] cutNode) {
    	System.out.print("sampled:\t[");
		for(int i=0; i < cutNode.length ; i++){
    		System.out.print(cutNode[i] + ",");
    	}	
    	System.out.println("]");

	}



	private double[] calculateConditionalProbabilityGivenLogNumeratorProb(
			double[] logNumeratorProb) {
		int numNodes = logNumeratorProb.length;
  		double maxLogProb = logNumeratorProb[0];
  		for(int i=0; i < numNodes; i++ ){
  			if(logNumeratorProb[i] > maxLogProb){
  				maxLogProb = logNumeratorProb[i];
  			}
  		}  		
  		
  		double sumLogDenominator = 0;
  		for(int i=0; i < numNodes; i++){
  			if(logNumeratorProb[i] != Double.NEGATIVE_INFINITY){
  				sumLogDenominator += Math.exp((logNumeratorProb[i]-maxLogProb));
  			}
  		}
  		sumLogDenominator = Math.log(sumLogDenominator) + maxLogProb;
  		
  		double sumProb = 0;
  		double []condProb = new double[numNodes]; 
  		for(int i=0; i < numNodes; i++){
  			condProb[i] = Math.exp( logNumeratorProb[i] - sumLogDenominator   );
			//System.out.println("condProb of site " + i + " = " + condProb[i]);
					sumProb +=condProb[i];
				if(condProb[i] > 0.01){
//					System.out.println("**site " + i + " with prob=" + condProb[i]  + "  steps from previous=" + numStepsFromOrigin[i]);
			}
  		}
  		return(condProb);
	}
	
	
	
	//Expect the ratio to be 1.
	private double verifyAssumption( double topLogMarginal,double[] logNumeratorProb) {
		
		
		for(int i=0; i < numNodes; i++){
			System.out.println(logNumeratorProb[i]);
		}
		int numNodes = logNumeratorProb.length;
  		double maxLogProb = logNumeratorProb[0];
  		for(int i=0; i < numNodes; i++ ){
  			if(logNumeratorProb[i] > maxLogProb){
  				maxLogProb = logNumeratorProb[i];
  			}
  		}  		
  		System.out.println("maxLogProb = " + maxLogProb);
  		
  		double sumLogDenominator = 0;
  		for(int i=0; i < numNodes; i++){
  			if(logNumeratorProb[i] != Double.NEGATIVE_INFINITY){
  				sumLogDenominator += Math.exp((logNumeratorProb[i]-maxLogProb));
  			}
  		}
  		System.out.println("tmp sum = " + sumLogDenominator);
  		sumLogDenominator = Math.log(sumLogDenominator) + maxLogProb;
  		
  		System.out.println("topLogMarginal = " + topLogMarginal);
  		System.out.println("sumLogDenominator = " + sumLogDenominator);
  		double ratio = Math.exp( Math.log(numNodes-1) + topLogMarginal - sumLogDenominator   );
  		
  		return(ratio);
	}




	private int checkSiteHasBeenAddedToOnIndicators(int curTest){
		int hasBeenAdded=0;
		for(int i=0; i < binSize; i++){
			if((int)status.getParameterValue(i) == 1 ){
				if( (int) breakPoints.getParameterValue(i) == curTest){
					hasBeenAdded=1;
					break;
				}
			}
		}
		return(hasBeenAdded);
	}


	private void setVirusLocationAndOffsets() {
		
		//change the mu in the toBin and fromBIn
		//borrow from getLogLikelihood:

		double[] meanYear = new double[binSize];
		double[] groupCount = new double[binSize];
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
					
		for(int i=0; i < binSize; i++){
			if(groupCount[i] > 0){
				meanYear[i] = meanYear[i]/groupCount[i];
			}
			//System.out.println(meanYear[i]);
		}


		mu0_offset = new double[binSize];
		//double[] mu1 = new double[maxLabel];
				
		
		//System.out.println("maxLabel=" + maxLabel);
		//now, change the mu..
		for(int i=0; i < binSize; i++){
			//System.out.println(meanYear[i]*beta);
			mu0_offset[i] =  meanYear[i];
			//System.out.println("group " + i + "\t" + mu0_offset[i]);
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
	   			//System.out.println("vloc="+ muValue + "," + muValue2);
    	}
    	
    	for(int i=0; i < numdata; i++){
        	int label = (int) clusterLabels.getParameterValue(i);
   			//if we want to apply the mean year virus cluster offset to the cluster
   			if(clusterOffsetsParameter != null){
   			//setting the clusterOffsets to be equal to the mean year of the virus cluster
   				// by doing this, the virus changes cluster AND updates the offset simultaneously
   				clusterOffsetsParameter.setParameterValue( i , mu0_offset[label]);
   			}
 				//		System.out.println("mu0_offset[label]=" + mu0_offset[label]);
 		//		System.out.println("clusterOffsets " +  i +" now becomes =" + clusterOffsetsParameter.getParameterValue(i) );   			
    	}

    	

    	
//    	System.out.println("===The on nodes===");
//    	for(int i=0; i < binSize; i++){	    
//    		if((int) excisionPoints.getParameterValue(i) == 1){
//    			System.out.println("Cluster node " + i + " = " + (int) indicators.getParameterValue(i) + "\tstatus=" + (int) excisionPoints.getParameterValue(i));
//    		}
//    	}
    	
		
	}




	private void setVirusLocationAutoCorrelatedModel() {
			int numNodes = treeModel.getNodeCount();
			double[][] nodeloc = new double[numNodes][2];
			
			//new - Trevor's autocorrelated model.
		//	System.out.println("Autocorrelated tree model");
			
			//given mu, excision points, and which ones are on...

			
			int[] nodeStatus = new int[numNodes];
			for(int i=0; i < numNodes; i ++){
				nodeStatus[i] = -1;
			}
			//convert to easy process format.
			for(int i=0; i < (binSize ); i++){
				if((int) status.getParameterValue(i) ==1){
					  nodeStatus[(int)breakPoints.getParameterValue(i)] = i;
				}
			}

			
//Testing:
//muValue[0] = 1;
//muValue2[0] = 1.5;	  
//nodeStatus[696] = 0;  
//muValue[1] = 10;
//muValue2[1] = 20;
//nodeStatus[607] = 1;  
//muValue[2] = 200;
//muValue2[2] = 300;

			
			//process the tree and get the vLoc of the viruses..
			//breadth first depth first..
			NodeRef cNode = treeModel.getRoot();
		    LinkedList<NodeRef> visitlist = new LinkedList<NodeRef>();
		    
		    visitlist.add(cNode);
		    
		    int countProcessed=0;
		    while(visitlist.size() > 0){
		    	
		    	
		    	countProcessed++;
		    	//assign value to the current node...
		    	if(treeModel.getParent(cNode) == null){
		    		Parameter curMu = mu.getParameter(0);
		    		nodeloc[cNode.getNumber()][0] =   curMu.getParameterValue(0);
		    		nodeloc[cNode.getNumber() ][1] = curMu.getParameterValue(1);
		    	}
		    	else{
		    		nodeloc[cNode.getNumber()][0] =   nodeloc[treeModel.getParent(cNode).getNumber()][0];
		    		nodeloc[cNode.getNumber()][1] =   nodeloc[treeModel.getParent(cNode).getNumber()][1];
		    		
		    		if(nodeStatus[cNode.getNumber()] != -1){
		    			//System.out.println("Run new location");
			    		Parameter curMu = mu.getParameter(nodeStatus[cNode.getNumber()] +1);
		    			nodeloc[cNode.getNumber()][0] += curMu.getParameterValue(0);
		    			nodeloc[cNode.getNumber()][1] += curMu.getParameterValue(1);	  			    			
		    		}
		    	}
		    	
		    	
		    	//add all the children to the queue
	  			for(int childNum=0; childNum < treeModel.getChildCount(cNode); childNum++){
	  				NodeRef node= treeModel.getChild(cNode,childNum);
	  				visitlist.add(node);
	  	        }
	  			
		  			
		  		visitlist.pop(); //now that we have finished visiting this node, pops it out of the queue
	
	  			if(visitlist.size() > 0){
	  				cNode = visitlist.getFirst(); //set the new first node in the queue to visit
	  			}
	  			
  			
		}
		    
		    //write the virus locations
		    for(int i=0; i < numdata; i++){
		    	Parameter vLocParameter = virusLocations.getParameter(i);
		    	vLocParameter.setParameterValue(0, nodeloc[membershipToClusterLabelIndexes[i]][0]);
		    	vLocParameter.setParameterValue(1, nodeloc[membershipToClusterLabelIndexes[i]][1]);
		    }
			
		    
		    //for(int i=0; i < numdata; i++){
				//Parameter vLocP= virusLocations.getParameter(i);
		    	//System.out.println("virus " + vLocP.getId() + "\t" + vLocP.getParameterValue(0) + "," + vLocP.getParameterValue(1)  );	  			    	
		    //}
		    
		    
		    
			//System.out.println("Processed " + countProcessed + " nodes");
			
			//System.out.println("Done");			
			
			//System.exit(0);		
	}





	private void relabelClusterLabels(int[] clusterLabel, int[] oldclusterLabel) {

    	int maxOldLabel = 0;
    	for(int i=0; i < oldclusterLabel.length; i++){
    		if(maxOldLabel < oldclusterLabel[i]){
    			maxOldLabel = oldclusterLabel[i];
    		}
    	}
    	
    	
        Map<Integer, Integer> m = new HashMap<Integer, Integer>();
        int[] isOldUsed = new int[ clusterLabel.length  ]; //an overkill - basically just need the max label in the old cluster
        
        for(int i=0; i < clusterLabel.length; i++){
        	
        	
    		if(m.get(new Integer(clusterLabel[i])) == null ){
    			if(isOldUsed[oldclusterLabel[i]] == 0){
    				m.put(new Integer(clusterLabel[i]), new Integer(oldclusterLabel[i]));
    				isOldUsed[oldclusterLabel[i]] = 1;
    			}
    			else{
    				maxOldLabel++;
    				m.put(new Integer(clusterLabel[i]), new Integer(maxOldLabel));
    			}
    			
    		}

    		clusterLabel[i] = m.get(new Integer( clusterLabel[i])).intValue();
    		
    	}
	}



	private int[] setClusterLabelsByTestCutNodeByNodeOrder(int[] testCutNode) {
        int []membership = determine_membershipByNodeOrder(treeModel, testCutNode, testCutNode.length);  // the time consuming step here.
        
  	   //The assumption that the first nodes being external node corresponding to the cluster labels IS FALSE
  	   //so I have to search for the matching indexes
       // for(int i=0; i < numdata; i++){
     	//   clusterLabels.setParameterValue( i, membership[membershipToClusterLabelIndexes[i]]);
       //}

        //to speed up the code
		int[] clusterLabel = new int[numdata];

        for(int i=0; i < numdata; i++){
        	clusterLabel[i] =  membership[membershipToClusterLabelIndexes[i]];
        }
        return(clusterLabel);
	}



	private void setMembershipToClusterLabelIndexes(){
        int numNodes = treeModel.getNodeCount();

  	   //I suspect this is an expensive operation, so I don't want to do it many times,
  	   //which is also unnecessary  - MAY have to update whenever a different tree is used.
         membershipToClusterLabelIndexes = new int[numdata]; 
         for(int i=0; i < numdata; i++){
  		   Parameter v = virusLocations.getParameter(i);
  		   String curName = v.getParameterName();
  		  // System.out.println(curName);
  		   int isFound = 0;
      	   for(int j=0; j < numNodes; j++){
      		   String treeId = treeModel.getTaxonId(j);
      		   if(curName.equals(treeId) ){
      		//	   System.out.println("  isFound at j=" + j);
      			   membershipToClusterLabelIndexes[i] = j;
      			   isFound=1;
      			   break;
      		   }	   
      	   }
      	   if(isFound ==0){
      		   System.out.println("not found. Exit now.");
      		   System.exit(0);
      	   }     	   
         }
    }
    
    //private void setClusterLabelsByTestCutNode(int[] testCutNode) {
    private int[] setClusterLabelsByTestCutNode(int[] testCutNode) {

        int []membership = determine_membership(treeModel, testCutNode, testCutNode.length);  // the time consuming step here.
        
  	   //The assumption that the first nodes being external node corresponding to the cluster labels IS FALSE
  	   //so I have to search for the matching indexes
       // for(int i=0; i < numdata; i++){
     	//   clusterLabels.setParameterValue( i, membership[membershipToClusterLabelIndexes[i]]);
       //}

        //to speed up the code
		int[] clusterLabel = new int[numdata];

        for(int i=0; i < numdata; i++){
        	clusterLabel[i] =  membership[membershipToClusterLabelIndexes[i]];
        }    
        
        return(clusterLabel);
		
	}



    
    //This function uses the breakPoints, but the doOperation doesn't use this anymore..
    //instead, it uses the cutNodes explicitly.
	private void setClusterLabels(int K_int) {

        int numNodes = treeModel.getNodeCount();
        int[] cutNodes = new int[K_int];
 	   int cutNum = 0;
 	   String content = "";
        for(int i=0; i < binSize; i++){
     	   if( (int) status.getParameterValue( i ) ==1 ){
     		   cutNodes[cutNum] = (int) breakPoints.getParameterValue(i);
     		   content += (int) breakPoints.getParameterValue(i) + ",";
     		   cutNum++;
     	   }
     	  
        }
       // System.out.println(content);
        
        if(cutNum != K_int){
        	System.out.println("cutNum != K_int. we got a problem");
        }
          
    //    for(int i=0; i < K_int; i++){
    // 	   System.out.println(cutNodes[i]);
     //   }
        
        //int []membership = determine_membership(treeModel, cutNodes, K_int-1);
        int []membership = determine_membership(treeModel, cutNodes, K_int);
        
        double uniqueCode = 0;
        for(int i=0; i < numNodes; i++){
        	uniqueCode += membership[i]*i;
        }
      //  System.out.println(" sum = " + uniqueCode);
        
     //   System.out.println("number of nodes = " + treeModel.getNodeCount());
      //  for(int i=0; i < treeModel.getNodeCount(); i++){
     //	   System.out.println(membership[i]);
      //  }
        
        
        //System.out.println("Done");
        
      //  for(int i=0; i < numdata; i++){
 	//	   Parameter v = virusLocations.getParameter(i);
 	//	   String curName = v.getParameterName();
 	//	   System.out.println("i=" + i + " = " + curName);       
 	//	}       
        
      //  for(int j=0; j < numdata; j++){
     //	   System.out.println("j=" + j + " = " + treeModel.getTaxonId(j));
      //  }
        
        
 	//   Parameter vv = virusLocations.getParameter(0);
 	 //  String curNamev = vv.getParameterName();
 	   
 	 //  System.out.println(curNamev + " and " +treeModel.getTaxonId(392) );
 	   //System.out.println(  curNamev.equals(treeModel.getTaxonId(392) )  );
 	   
        
        //System.exit(0);
        
 	  // System.out.println("numNodes=" + numNodes);
 	  // System.exit(0);
        //create dictionary:
 	   
 	   //I suspect this is an expensive operation, so I don't want to do it many times,
 	   //which is also unnecessary  - MAY have to update whenever a different tree is used.
        int []membershipToClusterLabelIndexes = new int[numdata]; 
        for(int i=0; i < numdata; i++){
 		   Parameter v = virusLocations.getParameter(i);
 		   String curName = v.getParameterName();
 		  // System.out.println(curName);
 		   int isFound = 0;
     	   for(int j=0; j < numNodes; j++){
     		   String treeId = treeModel.getTaxonId(j);
     		   if(curName.equals(treeId) ){
     		//	   System.out.println("  isFound at j=" + j);
     			   membershipToClusterLabelIndexes[i] = j;
     			   isFound=1;
     			   break;
     		   }
     		   
     	   }
     	   if(isFound ==0){
     		   System.out.println("not found. Exit now.");
     		   System.exit(0);
     	   }     	   
        }
        
        
       // System.exit(0);
        
      //  for(int i=0; i < numdata; i++){
     //	   System.out.println(membershipToClusterLabelIndexes[i]);
      //  }
       // System.exit(0);
        
        for(int i=0; i < numdata; i++){
     	   //The assumption that the first nodes being external node corresponding to the cluster labels IS FALSE
     	   //so I have to search for the matching indexes
     	   Parameter vloc = virusLocations.getParameter(i);
  
     	   
//must uncomment out because this sets the new partitioning ... now i am doing code testing.     	   
     	   clusterLabels.setParameterValue( i, membership[membershipToClusterLabelIndexes[i]]);
     	   //System.out.println(vloc.getParameterName() + " i="+ i + " membership=" + (int) clusterLabels.getParameterValue(i));
     	   
     	 //  Parameter v = virusLocations.getParameter(i);
     	  // System.out.println(v.getParameterName());
        }
        

    	
	}

    
    
    private static boolean isCutNode(int number, int cutNodes[], int numCut) {
    	if(numCut > 0){
    		for(int i=0; i < numCut; i++){
    			if(number == cutNodes[i]){
    				return true;
    			}
    		}
    	}
    	return false;
    }
    
    

  //traverse down the tree, top down, do calculation
  static int[] determine_membership(TreeModel treeModel, int[] cutNodes, int numCuts){

  	
  NodeRef root = treeModel.getRoot();

  int numClusters = 1;
  LinkedList<NodeRef> list = new LinkedList<NodeRef>();
  list.addFirst(root);

  int[] membership = new int[treeModel.getNodeCount()];
  for(int i=0; i < treeModel.getNodeCount(); i++){
  	membership[i] = -1;
  }
  membership[root.getNumber()] = 0; //root always given the first cluster
        
  while(!list.isEmpty()){
  	//do things with the current object
  	NodeRef curElement = list.pop();
  	//String content = "node #" + curElement.getNumber() +", taxon=" + treeModel.getNodeTaxon(curElement) + " and parent is = " ;
  	String content = "node #" + curElement.getNumber() +", taxon= " ;
  	if(treeModel.getNodeTaxon(curElement)== null){
  		content += "internal node\t";
  	}
  	else{
  		content += treeModel.getNodeTaxon(curElement).getId() + "\t";
  		//content += treeModel.getTaxonIndex(treeModel.getNodeTaxon(curElement)) + "\t";
  	}
  	
     	if(treeModel.getParent(curElement)== null){
  		//content += "no parent";
  	}
  	else{
  		//content += "parent node#=" + treeModel.getParent(curElement).getNumber();
  	}
  	
  	//cluster assignment:
  	if(!treeModel.isRoot(curElement)){
  	 if(isCutNode(curElement.getNumber(), cutNodes, numCuts)){
  	//if(isCutNode(curElement.getNumber())){
  		numClusters++ ;
  		membership[ curElement.getNumber() ] = numClusters - 1; 
    	}
  	else{
  		//inherit from parent's cluster assignment
  		membership[curElement.getNumber()] = membership[treeModel.getParent(curElement).getNumber()]; 
  	 }
  	        	
  	}//is not Root
  	content += " cluster = " + membership[curElement.getNumber()] ; 
  	
  //	System.out.println(content);

  	
      for(int childNum=0; childNum < treeModel.getChildCount(curElement); childNum++){
      	list.addFirst(treeModel.getChild(curElement,childNum));
      }
  }

   return(membership);
  }

    
  //traverse down the tree, top down, do calculation
  static int[] determine_membershipByNodeOrder(TreeModel treeModel, int[] cutNodes, int numCuts){


      Map<Integer, Integer> m = new HashMap<Integer, Integer>();
      for(int i=0; i < numCuts; i++){
    	  m.put(new Integer(cutNodes[i]), new Integer(i+1));
    	  
    	//  System.out.println(cutNodes[i] + "\t" + (i+1) );
      }
      
  	
  NodeRef root = treeModel.getRoot();

  int numClusters = 1;
  LinkedList<NodeRef> list = new LinkedList<NodeRef>();
  list.addFirst(root);

  int[] membership = new int[treeModel.getNodeCount()];
  for(int i=0; i < treeModel.getNodeCount(); i++){
  	membership[i] = -1;
  }
  membership[root.getNumber()] = 0; //root always given the first cluster
        
  while(!list.isEmpty()){
  	//do things with the current object
  	NodeRef curElement = list.pop();
  	//String content = "node #" + curElement.getNumber() +", taxon=" + treeModel.getNodeTaxon(curElement) + " and parent is = " ;
  	String content = "node #" + curElement.getNumber() +", taxon= " ;
  	if(treeModel.getNodeTaxon(curElement)== null){
  		content += "internal node\t";
  	}
  	else{
  		content += treeModel.getNodeTaxon(curElement).getId() + "\t";
  		//content += treeModel.getTaxonIndex(treeModel.getNodeTaxon(curElement)) + "\t";
  	}
  	
     	if(treeModel.getParent(curElement)== null){
  		//content += "no parent";
  	}
  	else{
  		//content += "parent node#=" + treeModel.getParent(curElement).getNumber();
  	}
  	
  	//cluster assignment:
  	if(!treeModel.isRoot(curElement)){
  	 if(isCutNode(curElement.getNumber(), cutNodes, numCuts)){
  	//if(isCutNode(curElement.getNumber())){
  		//numClusters++ ;
  		//membership[ curElement.getNumber() ] = numClusters - 1;
  		// System.out.println("get: curElement" + curElement.getNumber() + "\t" + m.get(new Integer( curElement.getNumber())));
  		membership[ curElement.getNumber()] = m.get(new Integer( curElement.getNumber()));
  		
    }
  	else{
  		//inherit from parent's cluster assignment
  		membership[curElement.getNumber()] = membership[treeModel.getParent(curElement).getNumber()]; 
  	 }
  	        	
  	}//is not Root
  	content += " cluster = " + membership[curElement.getNumber()] ; 
  	
  //	System.out.println(content);

  	
      for(int childNum=0; childNum < treeModel.getChildCount(curElement); childNum++){
      	list.addFirst(treeModel.getChild(curElement,childNum));
      }
  }

   return(membership);
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
		System.out.println("        	*      Rejected!");
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
    
	

           public final static String TREE_CLUSTERSEQUENTIAL_OPERATOR = "TreeClusterSequentialSampling";

              
            //MCMCOperator INTERFACE
            public final String getOperatorName() {
                return TREE_CLUSTERSEQUENTIAL_OPERATOR;
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
                
            	public final static String INDICATORS = "indicators";

                public final static String EXCISION_POINTS = "excisionPoints";


                public String getParserName() {
                    return TREE_CLUSTERSEQUENTIAL_OPERATOR;
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

                        cxo = xo.getChild(INDICATORS);
                        Parameter indicators = (Parameter) cxo.getChild(Parameter.class);
                      
                        cxo = xo.getChild(EXCISION_POINTS);
                        Parameter excisionPoints = (Parameter) cxo.getChild(Parameter.class);
                      
                    TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
 
                    AGLikelihoodTreeCluster agLikelihood = (AGLikelihoodTreeCluster) xo.getChild(AGLikelihoodTreeCluster.class);
                        
                    //return new ClusterAlgorithmOperator(virusLocations, mu, clusterLabels, k, weight, offsets, locationDrift, clusterOffsetsParameter);
                        return new TreeClusterSequentialSampling(virusLocations, mu, clusterLabels, k, weight, offsets,  clusterOffsetsParameter, indicators, excisionPoints, treeModel, agLikelihood);

                }

                //************************************************************************
                // AbstractXMLObjectParser implementation
                //************************************************************************

                public String getParserDescription() {
                    return "An operator that picks a new allocation of an item to a cluster under the Dirichlet process.";
                }

                public Class getReturnType() {
                    return TreeClusterSequentialSampling.class;
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
                       new ElementRule(INDICATORS, Parameter.class),
                       new ElementRule(EXCISION_POINTS, Parameter.class),
                       new ElementRule(TreeModel.class),

        
            };
            
            };


        
            public int getStepCount() {
                return 1;
            }

        }





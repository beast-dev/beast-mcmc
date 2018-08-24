package dr.evomodel.antigenic.phyloclustering.misc.obsolete;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.LinkedList;

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
public class TreeClusterGibbsOperator extends SimpleMCMCOperator implements GibbsOperator{

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
    
    private Parameter indicators = null;

    private int binSize=20;
 
    
    private Parameter excisionPoints;
    
    private TreeModel treeModel;

    
   // private int[] piIndicator = new int[numSites];

	
    //public ClusterAlgorithmOperator(MatrixParameter virusLocations, MatrixParameter mu, Parameter clusterLabels, Parameter K, double weight, Parameter virusOffsetsParameter, Parameter locationDrift_in, Parameter clusterOffsetsParameter) {
    public TreeClusterGibbsOperator(MatrixParameter virusLocations, MatrixParameter mu, Parameter clusterLabels, Parameter K, double weight, Parameter virusOffsetsParameter, Parameter clusterOffsetsParameter, Parameter indicatorsParameter, Parameter excisionPointsParameter, TreeModel treeModel_in, AGLikelihoodTreeCluster clusterLikelihood_in) {
    	
      	
    	System.out.println("Loading the constructor for ClusterGibbsOperator");
    	
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
    	this.indicators = indicatorsParameter;   
    	
    	this.excisionPoints = excisionPointsParameter;

        
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
    	
    	
    	int numNodes = treeModel.getNodeCount(); 

    	
    	double logHastingRatio = 0;  	
    	double chooseOperator = Math.random();
    	
    	int K_int = (int) K.getParameterValue(0);
    	    	
    	
    	int selectedI = -1;

    	
    	//System.out.println("AG likelihood cluster loaded");
    	//System.exit(0);
    	
    	//Gibbs move:
		double []logNumeratorProb = new double[numNodes];
        	

  		int isOn = 0;
  		int I_selected = -1;
  		
  		//obtain a random order of the "on" sites...
  		
  		while(isOn == 0){
  			I_selected = (int) (Math.floor(Math.random()*binSize));
  			isOn = (int) excisionPoints.getParameterValue(I_selected);
  	    	//find an "on" excision point to move.
  			if(isOn==1){
  				
//  				System.out.println("begin");
  				
  				int originalSite = (int) indicators.getParameterValue(I_selected);
  //				System.out.println("original site is = " + originalSite);
  				
  				
  				//Determining the number of steps from the original site
  				int []numStepsFromOrigin = new int[numNodes];
  				for(int i=0; i < numNodes; i++){
  					numStepsFromOrigin[i] = 100000;
  				}
  				
	  			int curElementNumber =(int) indicators.getParameterValue(I_selected);
	  			int rootElementNumber = curElementNumber;
	  			//System.out.println("curElementNumber=" + curElementNumber);
	  			NodeRef curElement = treeModel.getNode(curElementNumber); 
	  			
	  		    LinkedList<NodeRef> visitlist = new LinkedList<NodeRef>();
	  		    LinkedList<NodeRef> fromlist = new LinkedList<NodeRef>();
	  		    LinkedList<Integer> nodeLevel = new LinkedList<Integer>();
	  		    
	  		    LinkedList<Integer> possibilities = new LinkedList<Integer>();
	  		    
	  		    NodeRef dummyNode = null;
	  		    visitlist.add(curElement);
	  		    fromlist.add(dummyNode);
	  		    nodeLevel.add(new Integer(0));
	  		    
	  		    int maxNodeLevel = 1000;
  	  		    
  	  		  //System.out.println("root node " + curElement.getNumber());
  			    while(visitlist.size() > 0){
  					
  		  			if(treeModel.getParent(curElement) != null){
  		  				//add parent
  				  			NodeRef node= treeModel.getParent(curElement);	  		  			
  		  				if(fromlist.getFirst() != node){
  		  					if( nodeLevel.getFirst() < maxNodeLevel){
  		  						visitlist.add(node);
  		  		  				fromlist.add(curElement);
  		  		  				nodeLevel.add(new Integer(nodeLevel.getFirst()+1));
  		  		  				//System.out.println("node " +  node.getNumber() + " added, parent of " + curElement.getNumber());
  		  					}
  		  				}
  		  			}

  					
  		  			for(int childNum=0; childNum < treeModel.getChildCount(curElement); childNum++){
  		  				NodeRef node= treeModel.getChild(curElement,childNum);
  		  				if(fromlist.getFirst() != node){
  		  					if( nodeLevel.getFirst() < maxNodeLevel){
  		  						visitlist.add(node);
  		  						fromlist.add(curElement);
  		  						nodeLevel.add(new Integer(nodeLevel.getFirst()+1));
  		  						//System.out.println("node " +  node.getNumber() + " added, child of " + curElement.getNumber());
  		  					}
  		  				}
  		  	        }
  		  			



  		  			//System.out.println("visited " + curElement.getNumber());
  		  			//test if I can add curElement.getNumber()
  		  				int site_test = curElement.getNumber();
  			  			int hasBeenAdded=0;
  			  			for(int i=0; i < binSize; i++){
  			  				if( indicators.getParameterValue(i) == site_test){
  			  					hasBeenAdded=1;
  			  					break;
  			  				}
  			  			}
  			  			if(hasBeenAdded==0 || curElement.getNumber() == rootElementNumber ){
  			  				//System.out.println("to possibilities: add " + site_test);
  			  				numStepsFromOrigin[site_test] = nodeLevel.getFirst();
  			  				possibilities.addLast(new Integer( site_test));
  			  			}
  			  			else{
  			  				//System.out.println("element " + curElement.getNumber() + " is already an excision point");
  			  			}
  		  			
  	  		  			visitlist.pop();
  	  		  			fromlist.pop();
  	  		  			nodeLevel.pop();
  			  			
  		  			if(visitlist.size() > 0){
  		  				curElement = visitlist.getFirst();
  		  			}
  		  			
  		  			
  				}
  		  			
  				  				
  				
  				//Calculating the conditional probability
  				for(int curSite = 0; curSite < numNodes; curSite++){

					//check if a site has been added
					int hasBeenAdded=0;
		  			for(int i=0; i < binSize; i++){
		  				if( indicators.getParameterValue(i) == curSite){
		  					hasBeenAdded=1;
		  					break;
		  				}
		  			}
		  			//site that has been added will be zeroed out.
		  			if(hasBeenAdded==1){
		  				double inf = Double.NEGATIVE_INFINITY;
		  				logNumeratorProb[curSite]  = inf;

		  			}
		  			else{
		  				//calculate the numerator of the conditional probability

		  	    		//select that node, change the Clusterlabels, and virus offsets, calculate loglikelihood from AGLikelihoodCluster and store
		  	    		//perform Gibbs sampling
		  	    	   	//change the cluster labels and virus offsets
		  				int site_add = curSite; 
		  				//set to new sample
		  				indicators.setParameterValue(I_selected, site_add);


		  				
		  		    	//for each node that is not occupied,
		  				//select that node, change the Clusterlabels, and virus offsets, calculate loglikelihood from AGLikelihoodCluster and store
		  				//perform Gibbs sampling
		  			   	//change the cluster labels and virus offsets
		  			
		  				
		  		  		//swapped to the new on site
		  		  		//nonZeroIndexes[I_off] = site_on; 
		  		  		
		  		  		
		  		    	//K IS CHANGED ACCORDINGLY
		  				int K_count = 0; //K_int gets updated
		  				for(int i=0; i < binSize; i++){
		  					K_count += (int) excisionPoints.getParameterValue(i);
		  				}
		  				//System.out.println("K now becomes " + K_count);
		  				
		  				
		  				//Remove the commenting out later.
		  				K.setParameterValue(0, K_count); //update 
		  				K_int = K_count;
		  				
		  		    	//use the tree to re-partition according to the change.
		  				setClusterLabels(K_int);

		  					
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

		  	    	
		  			    	logNumeratorProb[curSite] = clusterLikelihood.getLogLikelihood();
		  			    	//System.out.println(clusterLikelihood.getLogLikelihood());
		  				
		  				
		  			}

  				}//close for loop
  				
  				selectedI = I_selected;
  			  	
  		  		
  		  		
  		  		double maxLogProb = logNumeratorProb[0];
  		  		for(int i=0; i < numNodes; i++ ){
  		  			if(logNumeratorProb[i] > maxLogProb){
  		  				maxLogProb = logNumeratorProb[i];
  		  			}
  		  		}
  		  		
  		  		//System.out.println(maxLogProb);
  		  		
  		  		
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
  //						System.out.println("**site " + i + " with prob=" + condProb[i]  + "  steps from previous=" + numStepsFromOrigin[i]);
  						
  					}
  		  		}
  		  		
  		  		//System.out.println("sum up to " + sumProb);
  		  		
  		  		
  		  		//int site_add = MathUtils.randomChoicePDF(condProb); //seems to not be working properly
  		  	int site_add = MathUtils.randomChoicePDF(condProb); //seems to not be working properly

				if(numStepsFromOrigin[site_add] >0){
					System.out.println("Gibbs move: indicator " + I_selected +" from site " + originalSite + " to "+  site_add + " , chosen with prob =" + condProb[site_add] + " steps from previous placement=" + numStepsFromOrigin[site_add] );
				}
				
//  		  		System.out.println("indicator " + I_selected +" from site " + originalSite + " to "+  site_add + " , chosen with prob =" + condProb[site_add] + " steps from previous placement=" + numStepsFromOrigin[site_add] );
  				indicators.setParameterValue(I_selected, site_add);
  				

  			} //close if isOn.
  		}    	  		
  		



 // 		MathUtils.randomChoicePDF(pairwiseDivergenceCountPerSite[group1][group2])
    	
  		
  	//	for(int curSite=0; curSite < numNodes; curSite++){
  	//		System.out.println("log condProb of site " + curSite + " = " + logNumeratorProb[curSite] + " and scaled=" +  (logNumeratorProb[curSite]-maxLogProb) );
  	//	}
  		
  		
		//normalize the score into probability distribution.
		
		//(int) Math.floor( Math.random()*numNodes );
//		int site_add = SAMPLE FROM THE DISTRIBUTION; 
		//set to new sample
	//	indicators.setParameterValue(I_selected, site_add);


		
    	//for each node that is not occupied,
		//select that node, change the Clusterlabels, and virus offsets, calculate loglikelihood from AGLikelihoodCluster and store
		//perform Gibbs sampling
	   	//change the cluster labels and virus offsets
	
 	
    	
		//====================================================================================================
		//After finishing the proposal
		//====================================================================================================
		
		
		
  		//swapped to the new on site
  		//nonZeroIndexes[I_off] = site_on; 
  		
  		
    	//K IS CHANGED ACCORDINGLY
		int K_count = 0; //K_int gets updated
		for(int i=0; i < binSize; i++){
			K_count += (int) excisionPoints.getParameterValue(i);
		}
		//System.out.println("K now becomes " + K_count);
		
		
		//Remove the commenting out later.
		K.setParameterValue(0, K_count); //update 
		K_int = K_count;
		
    	//use the tree to re-partition according to the change.
		setClusterLabels(K_int);

			
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

	    	

	    	
//	    	System.out.println("===The on nodes===");
//	    	for(int i=0; i < binSize; i++){	    
//	    		if((int) excisionPoints.getParameterValue(i) == 1){
//	    			System.out.println("Cluster node " + i + " = " + (int) indicators.getParameterValue(i) + "\tstatus=" + (int) excisionPoints.getParameterValue(i));
//	    		}
//	    	}
	    	
	    	
	    	
	    	
	    	//Hasting's Ratio is p(old |new)/ p(new|old)

	    	//System.out.println("Done doing operation!");
			
	    		    	
	    	
    	//return(logHastingRatio); //log hasting ratio
    	return(logHastingRatio);
    	
    }
    	
    	
    private void setClusterLabels(int K_int) {

        int numNodes = treeModel.getNodeCount();
        int[] cutNodes = new int[K_int];
 	   int cutNum = 0;
 	   String content = "";
        for(int i=0; i < binSize; i++){
     	   if( (int) excisionPoints.getParameterValue( i ) ==1 ){
     		   cutNodes[cutNum] = (int) indicators.getParameterValue(i);
     		   content += (int) indicators.getParameterValue(i) + ",";
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

  	
  	//TEMPORARY SOLUTION
      //load in the titer, corresponding to the taxon #.
  	
  	 TiterImporter titer = null ;
  	 
      FileReader fileReader;
  	try {
  		fileReader = new FileReader("/Users/charles/Documents/research/antigenic/GenoPheno/data/taxon_y_titer.txt");
  	     titer = new TiterImporter(fileReader);	

  	} catch (FileNotFoundException e) {
  		// TODO Auto-generated catch block
  		e.printStackTrace();
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
    
	

           public final static String TREE_CLUSTERGIBBS_OPERATOR = "TreeClusterGibbsOperator";

              
            //MCMCOperator INTERFACE
            public final String getOperatorName() {
                return TREE_CLUSTERGIBBS_OPERATOR;
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
                    return TREE_CLUSTERGIBBS_OPERATOR;
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
                        return new TreeClusterGibbsOperator(virusLocations, mu, clusterLabels, k, weight, offsets,  clusterOffsetsParameter, indicators, excisionPoints, treeModel, agLikelihood);

                }

                //************************************************************************
                // AbstractXMLObjectParser implementation
                //************************************************************************

                public String getParserDescription() {
                    return "An operator that picks a new allocation of an item to a cluster under the Dirichlet process.";
                }

                public Class getReturnType() {
                    return TreeClusterGibbsOperator.class;
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





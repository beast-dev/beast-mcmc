package dr.evomodel.antigenic.phyloclustering;

import java.util.LinkedList;

import dr.evolution.tree.NodeRef;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.MatrixParameter;
import dr.inference.model.Parameter;



public class TreeClusteringSharedRoutines {

	
	public static int[] setMembershipTreeToVirusIndexes(int numdata, MatrixParameter virusLocations, int numNodes, TreeModel treeModel ){

	  	   //I suspect this is an expensive operation, so I don't want to do it many times,
	  	   //which is also unnecessary  - MAY have to update whenever a different tree is used.
	         int[] correspondingTreeIndexForVirus = new int[numdata]; 
	         for(int i=0; i < numdata; i++){
	  		   Parameter v = virusLocations.getParameter(i);
	  		   String curName = v.getParameterName();
	  		  // System.out.println(curName);
	  		   int isFound = 0;
	      	   for(int j=0; j < numNodes; j++){
	      		   String treeId = treeModel.getTaxonId(j);
	      		   if(curName.equals(treeId) ){
	      		//	   System.out.println("  isFound at j=" + j);
	      			   correspondingTreeIndexForVirus[i] = j;
	      			   isFound=1;
	      			   break;
	      		   }	   
	      	   }
	      	   if(isFound ==0){
	      		   System.out.println("not found. Exit now.");
	      		   System.exit(0);
	      	   }     	   
	         }
	         
	         return(correspondingTreeIndexForVirus);
	    }
	
	public static void updateUndriftedVirusLocations(int numNodes, int numdata, TreeModel treeModel, MatrixParameter virusLocationsTreeNode, Parameter indicators, MatrixParameter mu, 		    MatrixParameter virusLocations, int[] correspondingTreeIndexForVirus){
			double[][] nodeloc = new double[numNodes][2];

			//process the tree and get the vLoc of the viruses..
			//breadth first depth first..
			NodeRef cNode = treeModel.getRoot();
		    LinkedList<NodeRef> visitlist = new LinkedList<NodeRef>();
		    
		    visitlist.add(cNode);
		    
		    int countProcessed=0;
		    while(visitlist.size() > 0){
		    	countProcessed++;
		    	//assign value to the current node...
		    	if(treeModel.getParent(cNode) == null){  //this means it is a root node
		    		Parameter curMu = mu.getParameter( cNode.getNumber() );
		    		//Parameter curMu = mu.getParameter(0);
		    		nodeloc[cNode.getNumber()][0]  = curMu.getParameterValue(0);
		    		nodeloc[cNode.getNumber() ][1] = curMu.getParameterValue(1);
		    		
		    		
		    		Parameter curVirusLoc = virusLocationsTreeNode.getParameter(cNode.getNumber());
		    		curVirusLoc.setParameterValue(0, curMu.getParameterValue(0) );
		    		curVirusLoc.setParameterValue(1, curMu.getParameterValue(1) );
		    	}
		    	else{
		    		nodeloc[cNode.getNumber()][0] =   nodeloc[treeModel.getParent(cNode).getNumber()][0];
		    		nodeloc[cNode.getNumber()][1] =   nodeloc[treeModel.getParent(cNode).getNumber()][1];
		    		
		    		if( (int) indicators.getParameterValue(cNode.getNumber()) == 1){
		    			Parameter curMu = mu.getParameter(cNode.getNumber() ); // no +1 because I don't need another mu- the root's mu takes care of the first cluster's mu 
			    		//Parameter curMu = mu.getParameter(cNode.getNumber() +1); //+1 because mu0 is reserved for the root.
		    			nodeloc[cNode.getNumber()][0] += curMu.getParameterValue(0);
		    			nodeloc[cNode.getNumber()][1] += curMu.getParameterValue(1);	  			    			
		    		}
		    		
		    		Parameter curVirusLoc = virusLocationsTreeNode.getParameter(cNode.getNumber());
		    		curVirusLoc.setParameterValue(0, nodeloc[cNode.getNumber()][0] );
		    		curVirusLoc.setParameterValue(1,nodeloc[cNode.getNumber()][1] );
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
		    	vLocParameter.setParameterValue(0, nodeloc[correspondingTreeIndexForVirus[i]][0]);
		    	vLocParameter.setParameterValue(1, nodeloc[correspondingTreeIndexForVirus[i]][1]);
		    }
			
		    
		    //for(int i=0; i < numdata; i++){
				//Parameter vLocP= virusLocations.getParameter(i);
		    	//System.out.println("virus " + vLocP.getId() + "\t" + vLocP.getParameterValue(0) + "," + vLocP.getParameterValue(1)  );	  			    	
		    //}
		    	
	}

	
	
	
	
	//may be very inefficient
	public static int findAnOnNodeIncludingRootRandomly(int numNodes, Parameter indicators) {
    	int isOn= 0;
    	int I_selected = -1;
  		while(isOn ==0){
  			I_selected = (int) (Math.floor(Math.random()*numNodes));
  			isOn = (int) indicators.getParameterValue(I_selected);  			
  		}    	  		
  		
    	return I_selected;
	}

	
	
	//Copied from TreeClusterAlgorithm  - should have put into the shared class...


	public static LinkedList<Integer> findActiveBreakpointsChildren(int selectedNodeNumber, int numNodes, TreeModel treeModel, Parameter indicators) {
			
			//a list of breakpoints...
			
			LinkedList<Integer> linkedList = new LinkedList<Integer>();
			int[] nodeBreakpointNumber = new int[numNodes];
						
			//int[] nodeStatus = new int[numNodes];
			//for(int i=0; i < numNodes; i ++){
			//	nodeStatus[i] = -1;
			//}
			
			//convert to easy process format.
			//for(int i=0; i < (binSize ); i++){
			//	if((int) indicators.getParameterValue(i) ==1){
			//		  nodeStatus[(int)breakPoints.getParameterValue(i)] = i;
			//	}
			//}
			
			//process the tree and get the vLoc of the viruses..
			//breadth first depth first..
			NodeRef cNode = treeModel.getRoot();
		    LinkedList<NodeRef> visitlist = new LinkedList<NodeRef>();

		    
		    visitlist.add(cNode);
		    
		    
		    //I am not sure if it still works......
		    
		    int countProcessed=0;
		    while(visitlist.size() > 0){
		    	
		    	
		    	countProcessed++;
		    	//assign value to the current node...
		    	if(treeModel.getParent(cNode) == null){
		    		//Parameter curMu = mu.getParameter(0);
		    		nodeBreakpointNumber[cNode.getNumber()] =   cNode.getNumber();
		    	}
		    	else{
		    		nodeBreakpointNumber[cNode.getNumber()] =   nodeBreakpointNumber[treeModel.getParent(cNode).getNumber()];
		    		//System.out.println("node#" + cNode.getNumber() + " is " + nodeBreakpointNumber[cNode.getNumber()]); 

		    		if( (int) indicators.getParameterValue(cNode.getNumber()) == 1){
		    			//System.out.println(cNode.getNumber() + " is a break point");
			    		//Parameter curMu = mu.getParameter(cNode.getNumber() +1); //+1 because mu0 is reserved for the root.
		    			//Parameter curMu = mu.getParameter(cNode.getNumber() ); //+1 because mu0 is reserved for the root.
			    		
			    		//see if parent's status is the same as the selectedIndex
			    		if( nodeBreakpointNumber[cNode.getNumber()] ==   selectedNodeNumber ){
			    			//System.out.println("hihi");
			    			linkedList.add( cNode.getNumber() );
			    		}
			    		//now, replace this nodeBreakpointNumber with its own node number
			    		nodeBreakpointNumber[cNode.getNumber()] = cNode.getNumber();
			    				    			  			    			
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
		    
		    //System.out.println("Now printing children of "  + selectedNodeNumber+":");
			//for(int i=0; i < linkedList.size(); i++){
			//	System.out.println( linkedList.get(i)  );
			//}
			
			return linkedList;
		}

		

	
}

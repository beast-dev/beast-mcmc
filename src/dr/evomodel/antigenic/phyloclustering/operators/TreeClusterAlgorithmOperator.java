package dr.evomodel.antigenic.phyloclustering.operators;


import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.Map;

import dr.evolution.tree.NodeRef;
import dr.evomodel.antigenic.AntigenicLikelihood;
import dr.evomodel.antigenic.phyloclustering.TreeClusteringSharedRoutines;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.MatrixParameter;
import dr.inference.model.Parameter;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;
import dr.math.distributions.MultivariateNormalDistribution;
import dr.util.DataTable;
import dr.xml.*;


/**
 * An operator to cluster viruses using a phylogenetic tree
 *
 * @author Charles Cheung
 * @author Trevor Bedford
 */
public class TreeClusterAlgorithmOperator extends SimpleMCMCOperator  {


    public final static String TREE_CLUSTERALGORITHM_OPERATOR = "treeClusterAlgorithmOperator";

	
	//Tuning parameters for proposals..
	private static final double WALK_SIZE = 4; // or 2 for +/- 1
    int maxNodeLevel = 4; //multistep - how many steps	
	
    
    
    //parameters
    private MatrixParameter mu = null;
    private Parameter clusterLabels = null;   
    private MatrixParameter virusLocations = null;
    private MatrixParameter serumLocations = null;

    private Parameter indicators;
    private Parameter muPrecision;
    private TreeModel treeModel;
  //  private AGLikelihoodTreeCluster clusterLikelihood = null;
    private AntigenicLikelihood clusterLikelihood = null; 


    private Parameter clusterLabelsTreeNode;
    private MatrixParameter virusLocationsTreeNode;
    
    private Parameter mu1Scale = null;
    private Parameter mu2Scale = null;
    private Parameter muMean = null;
   
    //-----------------------------------------------------------
    //I think these parameters are obsolete and should be removed 
	//private Parameter clusterOffsetsParameter;
	//private Parameter virusOffsetsParameter;
	//--------------------------------------------------
	
	private int numdata; //gets assigned in the constructor
	private int numNodes; //gets assigned in the constructor
	private int []correspondingTreeIndexForVirus = null; //relates treeModels's indexing system to cluster label's indexing system of viruses. Gets assigned

	//private int[] newClusterLabelArray; //for keeping the cluster labeling consistent
	//private int[] oldClusterLabelArray; //for keeping the cluster labeling consistent
	
	
    private int operatorSelect = -1; //keep track of which proposal gets called
     
    //For profiling acceptance rate
    private double []acceptNum;
    private double []rejectNum;
    
    private double []acceptDistance;
    private double []rejectDistance;
    
	private int moveCounter = 0; //counts how many moves have been proposed
	private int BURN_IN = 100000;
	private int frequencyPrintAcceptance = 1000000;
	
	//private int frequencyPrintActive = 10000; //for debugging, if printActiveNodes() is called.
	private int updateHotNodeFrequencey = 100000; //for the Propose_HotMultistepOnNodeFlipMu operator

	private double muDistance = -1; //for the Proposal_flipIBalanceRestrictive operator

	String[] operatorName =  {"Proposal_changeToAnotherNodeOn", 
							  "Proposal_changeMuFromPrior", 
							  "Proposal_flipIandChangeMu", 
							  "Proposal_changeAnOnMuWalk", 
							  "Proposal_multistepOnNode", 
							  "Propose_YandMu" , 
							  "Propose_YandI", 
							  "Propose_YandIandmu", 
							  "Propose_branchOffFlip", 
							  "Propose_multistepOnNodeFlipMu", 
							  "Propose_flipI", 
							  "Propose_changeOnMuAndBalance", 
							  "Proposal_changeMuFromWalk", 
							  "Proposal_changeAnOnMuFromPrior", 
							  "Propose_HotMultistepOnNodeFlipMu",
							  "Proposal_flipIBalance", 
							  "Proposal_OnMultistepIExchangeMuAndFlipAnotherI",  
							  "Proposal_changeRootMuWalk", 
							  "Proposal_changeRootMuWalkAndBalance",
							  "Proposal_flipIBalanceRestrictive"};
    
	
	
	//Decided after profiling acceptance..
	//Type:			highly crucial		Efficient	Booster to facilitate mixing
	//Exchange I:	0					9			14
	//Change mu:	1 (12)				11, 17
	//Flip I:		15								16
	//double[] operatorWeight = {1,1,0,0,0,0,0,0,0,1,0,1,0,0,1,1,1,0, 0.1};

	//Decided after profiling acceptance..
		//Type:			highly crucial		Efficient	Booster to facilitate mixing
		//Exchange I:	0					9			
		//Change mu:	1 (12)				11, 17
		//Flip I:		15								16
		double[] operatorWeight;// = {1,1,0,0,0,0,0,0,0,1,0,1,0,0,0,1,1,0, 0};
	//double[] operatorWeight = {0,1};
	
	
	//variables for the Propose_HotMultistepOnNodeFlipMu
    int[] hotNodes;
    int[] freqAcceptNode;
    private int curNode = 0;
    
    
    
    
    
    //Constructor
    public TreeClusterAlgorithmOperator(MatrixParameter virusLocations, 
    									MatrixParameter virusLocationsTreeNode_in,
    									MatrixParameter serumLocations, 
    									MatrixParameter mu, 
    									Parameter clusterLabels, 
    									double weight, 
    									//Parameter virusOffsetsParameter, 
    									//Parameter clusterOffsetsParameter, 
    									Parameter indicatorsParameter, 
    									TreeModel treeModel_in, 
    									AntigenicLikelihood clusterLikelihood_in,
    									Parameter muPrecision_in, 
    									DataTable<String[]> proposalWeightTable,
    									Parameter clusterLabelsTreeNode_in,
    									Parameter mu1Scale_in,
    									Parameter mu2Scale_in,
    									Parameter muMean_in) {
    	
    	operatorWeight = new double[proposalWeightTable.getRowCount()];
        for (int i = 0; i < proposalWeightTable.getRowCount(); i++) {
        	String[] values = proposalWeightTable.getRow(i);
        	operatorWeight[i] = Integer.parseInt(values[0]);
        }

    	
    	
    	acceptNum = new double[operatorWeight.length];
    	rejectNum = new double[operatorWeight.length];
    	for(int i=0; i < operatorWeight.length; i++){
    		acceptNum[i] = 0;
    		rejectNum[i] = 0;
    	}
    	
    	acceptDistance = new double[100];
    	rejectDistance = new double[100];
    	for(int i=0; i < 100; i++){
    		acceptDistance[i] = 0;
    		rejectDistance[i] = 0;
    	}
    	
    	
    	System.out.println("Loading the constructor for ClusterAlgorithmOperator");
		this.treeModel= treeModel_in;
    	this.mu = mu;
    	this.clusterLabels = clusterLabels;    	
        this.virusLocations = virusLocations;
        this.serumLocations = serumLocations;
       // this.virusOffsetsParameter = virusOffsetsParameter;
        //this.clusterOffsetsParameter = clusterOffsetsParameter;
    	this.indicators = indicatorsParameter;
    	this.clusterLikelihood = clusterLikelihood_in;
    	this.muPrecision = muPrecision_in;
    	this.clusterLabelsTreeNode = clusterLabelsTreeNode_in;
    	this.virusLocationsTreeNode = virusLocationsTreeNode_in;
    	
    	this.mu1Scale = mu1Scale_in;
    	this.mu2Scale = mu2Scale_in;
    	this.muMean = muMean_in;

    	
    	numNodes = treeModel.getNodeCount();
    	numdata = virusLocations.getColumnDimension();
    	//numdata = virusOffsetsParameter.getSize();
        System.out.println("numdata="+ numdata);
        
        setWeight(weight);
        
        System.out.println("Finished loading the constructor for ClusterAlgorithmOperator");
        
        
    	double sumOperatorWeight =0;
    	int numOp = operatorWeight.length;
    	for(int i=0; i < numOp; i++){
    		sumOperatorWeight += operatorWeight[i];
    	}
    	for(int i=0; i<numOp; i++){
    		operatorWeight[i] = operatorWeight[i]/sumOperatorWeight;
    	}
    	
    	System.out.println("#\tProposal\tCall Weight");
    	for(int i=0; i< numOp; i++){
    		System.out.println( i +"\t" + operatorName[i] + "\t" + operatorWeight[i]);
    	}
        
    	
    	//clusterLabelsTreeNode.setDimension(numNodes);
    	
    	correspondingTreeIndexForVirus = TreeClusteringSharedRoutines.setMembershipTreeToVirusIndexes(numdata, virusLocations, numNodes, treeModel);
    	//setMembershipTreeToVirusIndexes(); //run once to set up the 
    	//setClusterLabelsUsingIndicators(); 		//Update the cluster labels, after the indicators may have changed.
    	
    	TreeClusteringSharedRoutines.updateUndriftedVirusLocations(numNodes, numdata, treeModel, virusLocationsTreeNode, indicators, mu, virusLocations, correspondingTreeIndexForVirus);
    	//setVirusLocationAutoCorrelatedModel(); //set virus locations, given the breakpoints,status, and mu parameters	
    	//setClusterLabelsTreeNodesUsingIndicators();
    	CompositeSetClusterLabelsTreeNodesAndVirusesUsingIndicators();
		
		//to improve multistep sampling..
		hotNodes = new int[numNodes];
		freqAcceptNode = new int[numNodes];
    	for(int i=0; i < numNodes; i++){
    		hotNodes[i] = 1;
    		freqAcceptNode[i] = 0;
    	}
    
    	
  	  //NodeRef root = treeModel.getRoot();
  	  //System.out.println("Root node number = " + root.getNumber() );
  	  //System.exit(0);
    	
    	
//    	System.out.println("Hi");
 //   	for(int i=0; i < numNodes; i++){
  //  		System.out.println("node " + i + ": " + indicators.getParameterValue(i));
   // 	}
    	
    	
//    	CompositeSetClusterLabelsTreeNodesAndVirusesUsingIndicators();
    	
//    	for(int i=0; i < numNodes; i++){
 //  		   String treeId = treeModel.getTaxonId(i);
  //  		System.out.println("node " + i + " " +  treeId +": "+ clusterLabelsTreeNode.getParameterValue(i));
   // 	}
    	
   // 	System.out.println("=============================");
   // 	for(int i=0; i < numdata; i++){
 //		   	Parameter v = virusLocations.getParameter(i);
 //		   	String curName = v.getParameterName();
  //  		System.out.println( curName +": "+ clusterLabels.getParameterValue(i));
   // 	}
    	
    //	System.exit(0);
    	
    	
    //	loadClusterTreeNodes();
    //	setMembershipTreeToVirusIndexes(); //run once to set up the 
    //	PrintsetMembershipTreeToVirusIndexes();
    //	System.exit(0);
    	
    }
    


private void loadClusterTreeNodes() {

	FileReader fileReader2;
	try {
		//fileReader2 = new FileReader("/Users/charles/Documents/research/antigenic/GenoPheno/Gabriela/results/initialCondition/H3N2.serumLocs.log");
		//fileReader2 = new FileReader("/Users/charles/Documents/researchData/clustering/output/test25/run64/H3N2_mds.breakpoints.log");
	//	fileReader2 = new FileReader("/Users/charles/Documents/researchData/clustering/output/test25/run79/H3N2_mds.indicators.log");
		fileReader2 = new FileReader("/Users/charles/Documents/researchData/clustering/output/test26/run21/treeNodes120K.log");
		
	      BufferedReader bReader2 = new BufferedReader( fileReader2);

	      String line = null;

	      //skip to the last line
	      String testLine;
	      while ((testLine = bReader2.readLine()) != null){
	    	  line = testLine;
	      }

	    //  System.out.println(line);
	      
	      String datavalue[] = line.split("\t");

	      
	       //   System.out.println(serumLocationsParameter.getParameterCount());
	      for (int i = 0; i < treeModel.getNodeCount(); i++) {
	    	  
	    	  clusterLabelsTreeNode.setParameterValue(i, Double.parseDouble(datavalue[i+1]));
	    	 // System.out.println(datavalue[i*2+1]);
//	    	  System.out.println("indicator=" + indicators.getParameterValue(i));
	   	  
	      }
	      bReader2.close();
	
	} catch (FileNotFoundException e) {
		// TODO Auto-generated catch block
		e.printStackTrace();
	} catch (IOException e) {
		// TODO Auto-generated catch block
		e.printStackTrace();
	}        



}


    
    
 
    /**
	 * change the parameter and return the log hastings ratio.
     */
    public final double doOperation() {
    	double logHastingRatio = 0; //initiate the log Metropolis Hastings ratio of the MCMC
    	curNode = -1; //reset curNode. curNode is used to keep track of which node gets accepted... for the "hot" multistep proposal.
    
    	//Here, the tree doesn't change, so I don't need to repeat this procedure over and over again. just do it once in the constructor
    		//setMembershipTreeToVirusIndexes(); //run once in case the tree changes - to associate the tree with the virus indexes 
    		//numNodes = treeModel.getNodeCount(); 
  	    operatorSelect = MathUtils.randomChoicePDF(operatorWeight);
  	    
  	    
  	    // * * * * * * * * * * * * * * * * * *
  	    logHastingRatio = performProposal();  //This is the main routine for performing proposals.. it is broken down into many sub-routines to facilitate code maintenance
  	    // * * * * * * * * * * * * * * * * * *    

  	    //===  After the proposal, update cluster labels and virus locations ===
  	    // Note: some proposals may not involve the below steps.. For computational efficiency, I might want to only update if needed..  	    
  	   // setClusterLabelsUsingIndicators(); 		//1. Update the cluster labels, after the indicators parameters may have changed.
		
  	  TreeClusteringSharedRoutines.updateUndriftedVirusLocations(numNodes, numdata, treeModel, virusLocationsTreeNode, indicators, mu, virusLocations, correspondingTreeIndexForVirus);
  	    //setVirusLocationAutoCorrelatedModel(); //set virus locations, given the indicators and mu parameters
		//setClusterLabelsTreeNodesUsingIndicators();  
		CompositeSetClusterLabelsTreeNodesAndVirusesUsingIndicators();
		
		moveCounter ++;
			
    	return(logHastingRatio);    	
    }
    	



    private double performProposal() {
    	
		double logHR = 0;
		
		if(operatorSelect == 0){
			Proposal_changeToAnotherNodeOn();
   	 		}
 	    else if(operatorSelect == 1){
	 		logHR = Proposal_changeMuFromPrior(); //update mu 	   
 	    	}
 	    else if(operatorSelect == 2){
 	    	logHR = Proposal_flipIandChangeMu();
 	    	}
 	    else if(operatorSelect == 3){
 	    	logHR = Proposal_changeAnOnMuWalk();
 	    	}
 	    else if(operatorSelect == 4){
 	    	logHR = Proposal_multistepOnNode();
 	    	}
	 	else if(operatorSelect == 5){
		    logHR = Propose_YandMu();	  //NOT IMPLEMENTED COMPLETELY.. THE MH RATIO IS WRONG
	 		}
	 	else if(operatorSelect == 6){
	 		logHR = Propose_YandI();	 //NOT IMPLEMENTED COMPLETELY.. THE MH RATIO IS WRONG
		    }
	 	else if(operatorSelect == 7){
		    logHR = Propose_YandIandmu();	 //NOT IMPLEMENTED COMPLETELY.. THE MH RATIO IS WRONG
		    }
	 	else if(operatorSelect == 8){
		    logHR = Propose_branchOffFlip();	
		    }
	 	else if(operatorSelect == 9){
		    logHR = Propose_multistepOnNodeFlipMu();	
		    }
	 	else if(operatorSelect == 10){
	 		logHR = Proposal_flipI();	
		    }
	 	else if(operatorSelect == 11){
		    logHR = Propose_changeMuAndBalance();	
		    }
	 	else if(operatorSelect == 12){
	 		logHR = Proposal_changeMuWalk();
		    	
		    }
	 	else if(operatorSelect == 13){
		    logHR = Proposal_changeAnOnMuFromPrior();
		    }
	 	else if(operatorSelect == 14){
	 		logHR = Proposal_HotMultistepOnNodeFlipMu();
	 	}
	 	else if(operatorSelect == 15){
	 		logHR = Proposal_flipIBalance();
	 	}
	 	else if(operatorSelect == 16){
	 		logHR = Proposal_OnMultistepIExchangeMuAndFlipAnotherI(3);
	 	}
	 	else if(operatorSelect == 17){
	 		//System.out.println("hi: " + operatorWeight[18]);
	 		logHR = Proposal_changeRootMuWalk();
	 	}
	 	else if(operatorSelect == 18){
	 		logHR = Proposal_changeRootMuWalkAndBalance();
	 	}
		
	 	else if(operatorSelect == 19){
	 		logHR = Proposal_flipIBalanceRestrictive();
	 	}
		
	 	else if(operatorSelect == 100){
	 		test1();
		    }
	 	else if(operatorSelect == 101){
	 		test2();
		    }
	 	else if(operatorSelect == 102){	 		
	 		test3();
		    }
	 	else{
	 		//System.out.println("operatorSelect = " + operatorSelect);
	 		//System.out.println("Unimplemented operator. Quit now");
	 	}
		
		return(logHR);
	}


    
    
	//===============================================================================================
	//===============================================================================================
	
	//  BELOW IS A LIST OF PROPOSALS
	
	//===============================================================================================
	//===============================================================================================
	
	  
    
    
    
	private double Proposal_OnMultistepIExchangeMuAndFlipAnotherI(int maxNodeLevelHere) {
			
		double logHastingRatio = 0; 
		
		int rootNum = treeModel.getRoot().getNumber();
		
		//unlike the old version, self-move isn't allowed.
		
		int originalNode1 = findAnOnNodeRandomly();			//find an on-node
		
		
	//	System.out.print("Try " + originalNode1);
		
		int[] numStepsFromI_selected =determineTreeNeighborhood(originalNode1, 100000);
//System.out.print("[");		
		//1. Select an unoccupied site within some steps away from it.	 
		 LinkedList<Integer> possibilities1 = new LinkedList<Integer>();
		 for(int i=0; i < numNodes; i++){
			// System.out.println("#steps from I_selected " + numStepsFromI_selected[i]);
			 //make sure no self select
			 boolean isIn1 = numStepsFromI_selected[i] <= maxNodeLevelHere && numStepsFromI_selected[i] !=0 &&  i != rootNum;
			 if(isIn1){
				possibilities1.addLast(new Integer(i));
	//			System.out.print(i + ", ");
			 }
		 }//end for		
//System.out.println("]");		 
		
		 int numPossibilities1 = possibilities1.size();	
		 
		 if(numPossibilities1 > 0){
			 int whichMove = (int) (Math.floor(Math.random()*numPossibilities1)); //choose from possibilities
				 int site_add1 = possibilities1.get(whichMove).intValue();
			
				// System.out.println(" and select " + site_add1);
			   //  System.out.println("selected node = " + site_add1 + " that's " + numStepsFromI_selected[site_add1] + " steps from " + originalNode1);
				 
			indicators.setParameterValue(originalNode1, 0); // the existing indicator is now set to be off
			indicators.setParameterValue(site_add1, 1); //set the new selected index to the new node.
			
			curNode = site_add1;
			
		
			//Flip mu - so the neighbor that replaces the original node now also inherits the existing node's mu
			//Parameter originalNodeMu = mu.getParameter(originalNode1+1); //offset of 1
			Parameter originalNodeMu = mu.getParameter(originalNode1); 
			double[] tmp = originalNodeMu.getParameterValues();
		
			//Parameter newMu = mu.getParameter(site_add1+1); //offset of 1
			Parameter newMu = mu.getParameter(site_add1); 
			double[] tmpNew = newMu.getParameterValues();
			
			originalNodeMu.setParameterValue(0, tmpNew[0]);
			originalNodeMu.setParameterValue(1, tmpNew[1]);
		
			newMu.setParameterValue(0, tmp[0]);
			newMu.setParameterValue(1, tmp[1]);
		
			
			
		
			//new node calculation
			int[] numStepsNewNode =determineTreeNeighborhood(site_add1, 100000);
			 //System.out.print("[");
			//1. Select an unoccupied site within some steps away from it.
			 LinkedList<Integer> possibilities2 = new LinkedList<Integer>();
			 for(int i=0; i < numNodes; i++){
				// System.out.println("#steps from I_selected " + numStepsFromI_selected[i]);
				 //make sure no self select
				 boolean isIn2 = numStepsNewNode[i] <= maxNodeLevelHere && numStepsNewNode[i] !=0   && i != rootNum;
				 if(isIn2){
					possibilities2.addLast(new Integer(i));
				//	System.out.print(i + ", ");
				 }
	
			 }//end for
			 //System.out.println("]");
	
			 
			 int numPossibilities2 = possibilities2.size();		 
		
			
			
			//now need to combine the neighborhood of the pivot and the new pivot to determine the second move.
			LinkedList<Integer> jointPossibilities = new LinkedList<Integer>();
			int[] nodeInNeighborhood = new int[numNodes];
			 //System.out.print("[");
			for(int i=0; i < numPossibilities1; i++){
				int nodeNumber = possibilities1.get(i).intValue();
				jointPossibilities.addLast(new Integer(nodeNumber));
				nodeInNeighborhood[nodeNumber] = 1;
				//System.out.print(nodeNumber + ", ");
			}
	
			for(int i=0; i < numPossibilities2; i++){
				int nodeNumber = possibilities2.get(i).intValue();
				if(nodeInNeighborhood[nodeNumber] == 0){
					//add, since not in first list
					jointPossibilities.addLast(new Integer(nodeNumber));
				//	System.out.print(nodeNumber + ", ");
				}
			}
			 //System.out.println("]");
			
			int numJointPossibilities = jointPossibilities.size();
				
			//now, flip another multistep status:
			int whichMove2 = (int) (Math.floor(Math.random()*numJointPossibilities)); //choose from possibilities
			 int site_add2 = jointPossibilities.get(whichMove2).intValue();
		
			if((int)indicators.getParameterValue(site_add2) == 0 ){
				indicators.setParameterValue(site_add2, 1);
			}
			else{
				indicators.setParameterValue(site_add2, 0);
			}
				
	
			
		
		
		  //System.out.println("numPossibilities1=" + numPossibilities1 + " numPossibilities2 = " + numPossibilities2);

			//System.out.println("numJointPossibilities = " + numJointPossibilities);

	//forward: 1/N x 1/#from pivot  x 1/#from joint
	//backward: 1/N x 1/#from new pivot x 1/#from joint
	//backward/forward = 1/#from new pivot   / (1/#from pivot)
			logHastingRatio = Math.log( (1/ (double)numPossibilities2) / (1/ (double)numPossibilities1)  );
			
			//System.out.println("logHastingRatio = " + logHastingRatio);
			
		
			//System.out.println("need to test the code before using it");
			//System.exit(0);
		 }//if numPossibilities1 > 0
			return logHastingRatio;
	}



	private double Proposal_flipIBalance() {
		
		//System.out.println("hi it got run");
		//System.exit(0);
		
		

		//System.out.println("root: " + mu.getParameter(0).getParameterValue(0) + "," + mu.getParameter(0).getParameterValue(1));
		//for(int i=0; i < numNodes; i++){
			//if( (int) indicators.getParameterValue(i) == 1){
				//System.out.println(i + ": " + mu.getParameter(i+1).getParameterValue(0) + "," + mu.getParameter(i+1).getParameterValue(1));
			//}
		//}
		
		
		int node = findNodeRandomly();
		//int node = (int) (Math.floor(Math.random()*numNodes));
	//node = 785;
		//System.out.println("selected node " + node);
		//double[] originalValues = mu.getParameter(node +1).getParameterValues();	
		double[] originalValues = mu.getParameter(node ).getParameterValues();
		//System.out.println(originalValues[0] + " and " + originalValues[1]);
		//a. by turning on the selected node, each child of this node should be updated to keep the absolute location of 
		//the child cluster fixed as before
		LinkedList<Integer> childrenOriginalNode = findActiveBreakpointsChildren(node);		
		
		if((int)indicators.getParameterValue(node) == 0 ){
			indicators.setParameterValue(node, 1);
			//System.out.println("turn it on");

			for(int i=0; i < childrenOriginalNode.size(); i++){
				int muIndexNum = childrenOriginalNode.get(i).intValue() ;
				//int muIndexNum = childrenOriginalNode.get(i).intValue() + 1;
				Parameter curMu = mu.getParameter( muIndexNum );
				double curMu_original0 = curMu.getParameterValue( 0);
				mu.getParameter(muIndexNum).setParameterValue(0, curMu_original0 - originalValues[0]);
				double curMu_original1 = curMu.getParameterValue( 1);
				mu.getParameter(muIndexNum).setParameterValue(1, curMu_original1 - originalValues[1]);
				//System.out.println( " " + ( muIndexNum - 1) + " is a child");
			}

		}
		else{
			indicators.setParameterValue(node, 0);
			//System.out.println("turn it off");
			for(int i=0; i < childrenOriginalNode.size(); i++){
				int muIndexNum = childrenOriginalNode.get(i).intValue() ;
				//int muIndexNum = childrenOriginalNode.get(i).intValue() + 1;
				Parameter curMu = mu.getParameter( muIndexNum );
				double curMu_original0 = curMu.getParameterValue( 0);
				mu.getParameter(muIndexNum).setParameterValue(0, curMu_original0 + originalValues[0]);
				double curMu_original1 = curMu.getParameterValue( 1);
				mu.getParameter(muIndexNum).setParameterValue(1, curMu_original1 + originalValues[1]);
				//System.out.println( " " + ( muIndexNum - 1) + " is a child");
			}
			
		}
		
		double coord1 = mu1Scale.getParameterValue(0)*originalValues[0];
		double coord2 = mu2Scale.getParameterValue(0)*originalValues[1];
	
		muDistance = Math.sqrt( coord1*coord1 + coord2*coord2);
		//System.out.println("root: " + mu.getParameter(0).getParameterValue(0) + "," + mu.getParameter(0).getParameterValue(1));
		//for(int i=0; i < numNodes; i++){
			//if( (int) indicators.getParameterValue(i) == 1){
				//System.out.println(i + ": " + mu.getParameter(i+1).getParameterValue(0) + "," + mu.getParameter(i+1).getParameterValue(1));
			//}
		//}
		
		
		//System.exit(0);
		
		return(0);
	}

	
	private double Proposal_flipIBalanceRestrictive(){
		int node = findRestrictedNodeRandomly(2); //neighborhood size is 2
		if(node != -1){
			double[] originalValues = mu.getParameter(node ).getParameterValues();
			//System.out.println(originalValues[0] + " and " + originalValues[1]);
			//a. by turning on the selected node, each child of this node should be updated to keep the absolute location of 
			//the child cluster fixed as before
			LinkedList<Integer> childrenOriginalNode = findActiveBreakpointsChildren(node);		
			
			if((int)indicators.getParameterValue(node) == 0 ){
				indicators.setParameterValue(node, 1);
				for(int i=0; i < childrenOriginalNode.size(); i++){
					int muIndexNum = childrenOriginalNode.get(i).intValue() ;
					//int muIndexNum = childrenOriginalNode.get(i).intValue() + 1;
					Parameter curMu = mu.getParameter( muIndexNum );
					double curMu_original0 = curMu.getParameterValue( 0);
					mu.getParameter(muIndexNum).setParameterValue(0, curMu_original0 - originalValues[0]);
					double curMu_original1 = curMu.getParameterValue( 1);
					mu.getParameter(muIndexNum).setParameterValue(1, curMu_original1 - originalValues[1]);
					//System.out.println( " " + ( muIndexNum - 1) + " is a child");
				}
	
			}
			else{
				indicators.setParameterValue(node, 0);
				//System.out.println("turn it off");
				for(int i=0; i < childrenOriginalNode.size(); i++){
					int muIndexNum = childrenOriginalNode.get(i).intValue() ;
					//int muIndexNum = childrenOriginalNode.get(i).intValue() + 1;
					Parameter curMu = mu.getParameter( muIndexNum );
					double curMu_original0 = curMu.getParameterValue( 0);
					mu.getParameter(muIndexNum).setParameterValue(0, curMu_original0 + originalValues[0]);
					double curMu_original1 = curMu.getParameterValue( 1);
					mu.getParameter(muIndexNum).setParameterValue(1, curMu_original1 + originalValues[1]);
					//System.out.println( " " + ( muIndexNum - 1) + " is a child");
				}
				
			}			
		
			muDistance = Math.sqrt( originalValues[0]*originalValues[0] + originalValues[1]*originalValues[1]);
			return(0);
		}
		else{
			//don't accept the move, since no valid choice
			return(Double.NEGATIVE_INFINITY);
		}
	}
	
	private int findRestrictedNodeRandomly(double neighborhood) {
		
		double mu1ScaleValue = mu1Scale.getParameterValue(0);
		double mu2ScaleValue = mu2Scale.getParameterValue(0);
		

		int rootNode = treeModel.getRoot().getNumber();

		int numQualified = 0;
		int[] qualifiedNodes = new int[numNodes];
		for(int i=0; i < numNodes; i++){
			if(i != rootNode){
				Parameter curMu = mu.getParameter(i);
				double coord1 = mu1ScaleValue * curMu.getParameterValue(0); 
				double coord2 = mu2ScaleValue * curMu.getParameterValue(1);
				
		//		double coord1 = curMu.getParameterValue(0);
			//	double coord2 = curMu.getParameterValue(1);
				double dist = Math.sqrt( coord1*coord1 + coord2*coord2);
				if(dist < neighborhood){
					qualifiedNodes[numQualified] = i;
					numQualified++;
				}
			}
		}

		//now draw 
		if( numQualified >0){
			int ranSelect = (int) (Math.random()*numQualified); 	
			int selectedNode = qualifiedNodes[ranSelect];
			return selectedNode;
		}
		return -1; // no node qualified, return -1
			
	}



	private double Proposal_changeRootMuWalkAndBalance(){
		
		int rootNum = treeModel.getRoot().getNumber();
		
		int dimSelect = (int) Math.floor( Math.random()* 2 );   		  	    		
		double change = Math.random()*WALK_SIZE- WALK_SIZE/2 ; 	
		//double originalValue = mu.getParameter(0).getParameterValue(dimSelect);		
		//mu.getParameter(0).setParameterValue(dimSelect, originalValue + change);
		double originalValue = mu.getParameter(rootNum).getParameterValue(dimSelect);
		mu.getParameter(rootNum).setParameterValue(dimSelect, originalValue + change);
		
		//a. by removing the selected node, each child of this node should be updated to keep the absolute location of 
		//the child cluster fixed as before
		//LinkedList<Integer> childrenOriginalNode = findActiveBreakpointsChildren(-1); //find the root's children		
		LinkedList<Integer> childrenOriginalNode = findActiveBreakpointsChildren(rootNum); //find the root's children
		for(int i=0; i < childrenOriginalNode.size(); i++){
			//int muIndexNum = childrenOriginalNode.get(i).intValue() + 1;
			int muIndexNum = childrenOriginalNode.get(i).intValue() ;
			Parameter curMu = mu.getParameter( muIndexNum );
			double curMu_original = curMu.getParameterValue( dimSelect);
			mu.getParameter(muIndexNum).setParameterValue(dimSelect, curMu_original - change);
			
			//System.out.println( " " + ( muIndexNum - 1) + " is a child");
			
		}

		
		return(0);
	
	}


	private double Propose_changeMuAndBalance() {

		//System.out.println("root: " + mu.getParameter(0).getParameterValue(0) + "," + mu.getParameter(0).getParameterValue(1));
		//for(int i=0; i < numNodes; i++){
		//	if( (int) indicators.getParameterValue(i) == 1){
		//		System.out.println(i + ": " + mu.getParameter(i+1).getParameterValue(0) + "," + mu.getParameter(i+1).getParameterValue(1));
		//	}
		//}
		
		//first, randomly select an "on" node to overwrite
		int originalNode = findAnOnNodeIncludingRootRandomly();			//find an on-node	
		//originalNode = 673;
		
		
		//if(originalNode == 802){
			//System.out.println(treeModel.getRoot().getNumber());
			//System.out.println("I am walking 802!");
		//}
		//unbounded walk
		int dimSelect = (int) Math.floor( Math.random()* 2 );   		  	    		
		double change = Math.random()*WALK_SIZE- WALK_SIZE/2 ; 	
		
		
		//dimSelect = 0;
		//change = 10;	
		//double originalValue = mu.getParameter(originalNode +1).getParameterValue(dimSelect);
		double originalValue = mu.getParameter(originalNode).getParameterValue(dimSelect);
		//System.out.println("originalValue = " + originalValue);
		mu.getParameter(originalNode ).setParameterValue(dimSelect, originalValue + change);
		//mu.getParameter(originalNode + 1).setParameterValue(dimSelect, originalValue + change);
  		//System.out.println("original node = " + originalNode);
	
		//a. by removing the selected node, each child of this node should be updated to keep the absolute location of 
		//the child cluster fixed as before
		LinkedList<Integer> childrenOriginalNode = findActiveBreakpointsChildren(originalNode);	
		//if(originalNode == 802){
			//System.out.println("number of child = " + childrenOriginalNode.size());
		//}
		for(int i=0; i < childrenOriginalNode.size(); i++){
			//int muIndexNum = childrenOriginalNode.get(i).intValue() + 1;
			int muIndexNum = childrenOriginalNode.get(i).intValue() ;
			//if(originalNode == 802){
				//System.out.println(" " + muIndexNum + " is a child");
			//}
			Parameter curMu = mu.getParameter( muIndexNum );
			double curMu_original = curMu.getParameterValue( dimSelect);
			mu.getParameter(muIndexNum).setParameterValue(dimSelect, curMu_original - change);
			
			//System.out.println( " " + ( muIndexNum - 1) + " is a child");
			
		}
		
	//	System.out.println("root: " + mu.getParameter(0).getParameterValue(0) + "," + mu.getParameter(0).getParameterValue(1));
	//	for(int i=0; i < numNodes; i++){
	//		if( (int) indicators.getParameterValue(i) == 1){
	//			System.out.println(i + ": " + mu.getParameter(i+1).getParameterValue(0) + "," + mu.getParameter(i+1).getParameterValue(1));
	//		}
	//	}
		
		
		return(0);
	}



	private double Proposal_HotMultistepOnNodeFlipMu() {
		
		int rootNum = treeModel.getRoot().getNumber();

		//unlike the old version, self-move isn't allowed.
		
		int originalNode1 = findAnOnNodeRandomly();			//find an on-node
		
		//System.out.print("Try " + originalNode1);
		int[] numStepsFromI_selected =determineTreeNeighborhood(originalNode1, 100000);
		
		//1. Select an unoccupied site within some steps away from it.	 
		 LinkedList<Integer> possibilities1 = new LinkedList<Integer>();
		 for(int i=0; i < numNodes; i++){
			// System.out.println("#steps from I_selected " + numStepsFromI_selected[i]);
			 //make sure no self select
			 boolean isIn1 = numStepsFromI_selected[i] <= maxNodeLevel && numStepsFromI_selected[i] !=0  &&  i != rootNum;
			 if(isIn1 && hotNodes[i] ==1){
				possibilities1.addLast(new Integer(i));
			 }
		 }//end for		
		 
		 int numPossibilities1 = possibilities1.size();		
		 
		 //if there is a single legal configuration to switch to
		 if(numPossibilities1 > 0){
		 int whichMove = (int) (Math.floor(Math.random()*numPossibilities1)); //choose from possibilities
			 int site_add1 = possibilities1.get(whichMove).intValue();
		//	 System.out.println(" and select " + site_add1);
		// System.out.println("selected node = " + site_add1 + " that's " + numStepsFromI_selected[site_add1] + " steps from " + originalNode1);
			 
		indicators.setParameterValue(originalNode1, 0); // the existing indicator is now set to be off
		indicators.setParameterValue(site_add1, 1); //set the new selected index to the new node.
		
		

		//Flip mu - so the neighbor that replaces the original node now also inherits the existing node's mu
		//Parameter originalNodeMu = mu.getParameter(originalNode1+1); //offset of 1
		Parameter originalNodeMu = mu.getParameter(originalNode1); 
		double[] tmp = originalNodeMu.getParameterValues();

		//Parameter newMu = mu.getParameter(site_add1+1); //offset of 1
		Parameter newMu = mu.getParameter(site_add1); 
		double[] tmpNew = newMu.getParameterValues();
		
		originalNodeMu.setParameterValue(0, tmpNew[0]);
		originalNodeMu.setParameterValue(1, tmpNew[1]);

		newMu.setParameterValue(0, tmp[0]);
		newMu.setParameterValue(1, tmp[1]);
		
	
		//backward calculation
		int[] numStepsBackward =determineTreeNeighborhood(site_add1, 100000);
		
		//1. Select an unoccupied site within some steps away from it.
		 LinkedList<Integer> possibilities2 = new LinkedList<Integer>();
		 for(int i=0; i < numNodes; i++){
			// System.out.println("#steps from I_selected " + numStepsFromI_selected[i]);
			 //make sure no self select
			 boolean isIn2 = numStepsBackward[i] <= maxNodeLevel && numStepsBackward[i] !=0   &&  i != rootNum;
			 if(isIn2 && hotNodes[i] ==1){
				possibilities2.addLast(new Integer(i));
			 }
		 }//end for
		 int numPossibilities2 = possibilities2.size();		 
		 //	 System.out.println("numPossibilities1=" + numPossibilities1 + " numPossibilities2 = " + numPossibilities2);

		 
  		double logHastingRatio = Math.log( (1/ (double)numPossibilities2) / (1/ (double)numPossibilities1)  );
  		//System.out.println("logHastingRatio = " + logHastingRatio);
		return logHastingRatio;
		 }
		 else{
			 return Double.NEGATIVE_INFINITY;
		 }
	}


	
	private double Propose_multistepOnNodeFlipMu() {
		double logHastingRatio = 0;
		
		//unlike the old version, self-move isn't allowed.
		int rootNum = treeModel.getRoot().getNumber();
		
		int originalNode1 = findAnOnNodeRandomly();			//find an on-node
		
		//System.out.print("Try " + originalNode1);
		int[] numStepsFromI_selected =determineTreeNeighborhood(originalNode1, 100000);
		
		//1. Select an unoccupied site within some steps away from it.	 
		 LinkedList<Integer> possibilities1 = new LinkedList<Integer>();
		 for(int i=0; i < numNodes; i++){
			// System.out.println("#steps from I_selected " + numStepsFromI_selected[i]);
			 //make sure no self select
			 boolean isIn1 = numStepsFromI_selected[i] <= maxNodeLevel && numStepsFromI_selected[i] !=0 &&  i != rootNum
					 && (int) indicators.getParameterValue(i) == 0;
			 if(isIn1){
				possibilities1.addLast(new Integer(i));
			 }
		 }//end for		
		 
		 int numPossibilities1 = possibilities1.size();	
		 if(numPossibilities1 > 0){
			 int whichMove = (int) (Math.floor(Math.random()*numPossibilities1)); //choose from possibilities
				 int site_add1 = possibilities1.get(whichMove).intValue();
			//	 System.out.println(" and select " + site_add1);
			// System.out.println("selected node = " + site_add1 + " that's " + numStepsFromI_selected[site_add1] + " steps from " + originalNode1);
				 
			indicators.setParameterValue(originalNode1, 0); // the existing indicator is now set to be off
			indicators.setParameterValue(site_add1, 1); //set the new selected index to the new node.
			
			curNode = site_add1;
			
	
			//Flip mu - so the neighbor that replaces the original node now also inherits the existing node's mu
			//Parameter originalNodeMu = mu.getParameter(originalNode1+1); //offset of 1
			Parameter originalNodeMu = mu.getParameter(originalNode1);
			double[] tmp = originalNodeMu.getParameterValues();
	
			//Parameter newMu = mu.getParameter(site_add1+1); //offset of 1
			Parameter newMu = mu.getParameter(site_add1); 
			double[] tmpNew = newMu.getParameterValues();
			
			originalNodeMu.setParameterValue(0, tmpNew[0]);
			originalNodeMu.setParameterValue(1, tmpNew[1]);
	
			newMu.setParameterValue(0, tmp[0]);
			newMu.setParameterValue(1, tmp[1]);
			
		
			//backward calculation
			int[] numStepsBackward =determineTreeNeighborhood(site_add1, 100000);
			
			//1. Select an unoccupied site within some steps away from it.
			 LinkedList<Integer> possibilities2 = new LinkedList<Integer>();
			 for(int i=0; i < numNodes; i++){
				// System.out.println("#steps from I_selected " + numStepsFromI_selected[i]);
				 //make sure no self select
				 boolean isIn2 = numStepsBackward[i] <= maxNodeLevel && numStepsBackward[i] !=0   &&  i != rootNum
						 && (int) indicators.getParameterValue(i) == 0;
				 if(isIn2){
					possibilities2.addLast(new Integer(i));
				 }
			 }//end for
			 int numPossibilities2 = possibilities2.size();		 
		 //	 System.out.println("numPossibilities1=" + numPossibilities1 + " numPossibilities2 = " + numPossibilities2);
		 	if(numPossibilities2 > 0){
		 		logHastingRatio = Math.log( (1/ (double)numPossibilities2) / (1/ (double)numPossibilities1)  );
		 	}
		 }
  		//System.out.println("logHastingRatio = " + logHastingRatio);
		return logHastingRatio;
	}



	private double Propose_branchOffFlip() {
		
		
		//first, randomly select an "on" node to overwrite
		int originalNode = findAnOnNodeRandomly();			//find an on-node	
			
			
		//second, randomly select a destination
		int site_add = findAnOffNodeRandomly();			//sample a node that's not in the cluster.
	
		//existing mu
		Parameter selectedMu = mu.getParameter(originalNode ) ;
		//Parameter selectedMu = mu.getParameter(originalNode +1) ;
			
		double selectedMu0 = selectedMu.getParameterValue(0);
		double selectedMu1 = selectedMu.getParameterValue(1);
		
			
		//a. by removing the selected node, each child of this node should be updated to keep the absolute location of 
		//the child cluster fixed as before
		LinkedList<Integer> childrenOriginalNode = findActiveBreakpointsChildren(originalNode);	
		
		for(int i=0; i < childrenOriginalNode.size(); i++){
			int muIndexNum = childrenOriginalNode.get(i).intValue() ;
			//int muIndexNum = childrenOriginalNode.get(i).intValue() + 1;
			Parameter curMu = mu.getParameter( muIndexNum );
			double mu0 = curMu.getParameterValue(0) + selectedMu0;
			double mu1 = curMu.getParameterValue(1) + selectedMu1;
			mu.getParameter(muIndexNum).setParameterValue(0, mu0);
			mu.getParameter(muIndexNum).setParameterValue(1, mu1);

		}


			//set indicators AND NEW MU
			indicators.setParameterValue(site_add, 1);
			//indicators.setParameterValue(originalNode, 0); //just flip it on. do not replace

			//I think this generate a situation where if a mu walks off, then it gets the breakpoint that doesn't partition.
			//this creates a scenario where a breakpoint is lost
			//double change = Math.random()*WALK_SIZE - WALK_SIZE ; 
			//double newMu0 = selectedMu0 + change;
			//double change2 = Math.random()*WALK_SIZE- WALK_SIZE ; 
			//double newMu1 = selectedMu1 + change2;
			
			

			double[] oldValues = mu.getParameter(site_add).getParameterValues();
			//double[] oldValues = mu.getParameter(site_add+1).getParameterValues();
			
			//System.out.println(oldValues[0]  + ", " + oldValues[1]);
			//instead, sample from the normal distribution
			double[] mean = new double[2];
			mean[0] = 0;
			mean[1] = 0;
			double[][] precisionM = new double[2][2];
			//double precision = 1/TreeClusterViruses.getSigmaSq();
			double precision = muPrecision.getParameterValue(0);
			precisionM[0][0] = precision;
			precisionM[0][1] = 0;
			precisionM[1][0] = 0;
			precisionM[1][1] = precision;
			
			
			double[] values = MultivariateNormalDistribution.nextMultivariateNormalPrecision(mean, precisionM);
			//System.out.println(values[0]  + ", " + values[1]); 
			
			mu.getParameter(site_add).setParameterValue(0,values[0]);
			mu.getParameter(site_add).setParameterValue(1,values[1]);

			//mu.getParameter(site_add+1).setParameterValue(0,values[0]);
			//mu.getParameter(site_add+1).setParameterValue(1,values[1]);

			

			selectedMu.setParameterValue(0, values[0]);
			selectedMu.setParameterValue(1, values[1]);
			
			
			

			//b. by adding the new selected node, each child of this new node should be updated to keep the absolute location of 
		//the child cluster fixed as before
			LinkedList<Integer> childrenNewNode = findActiveBreakpointsChildren(site_add);
			
			for(int i=0; i < childrenNewNode.size(); i++){
				//int muIndexNum = childrenNewNode.get(i).intValue() + 1;
				int muIndexNum = childrenNewNode.get(i).intValue() ;
				Parameter curMu = mu.getParameter( muIndexNum);
				double mu0 = curMu.getParameterValue(0) - values[0];
				double mu1 = curMu.getParameterValue(1) - values[1];
				mu.getParameter(muIndexNum).setParameterValue(0, mu0);
				mu.getParameter(muIndexNum).setParameterValue(1, mu1);
			}
			
			
			
			
			
			double logHastingRatio = MultivariateNormalDistribution.logPdf(oldValues, mean, precision, 1) - MultivariateNormalDistribution.logPdf(values, mean, precision, 1) ; 


			return(logHastingRatio);
	}



	private double Propose_YandIandmu() {

		int rootNum = treeModel.getRoot().getNumber();

		//first, find a random Y and walk
			int serum_selected = (int) (Math.floor(Math.random()*getNumSera()));
			
			MatrixParameter serumLocations = getSerumLocationsParameter();
			Parameter serum = serumLocations.getParameter(serum_selected);
			int whichDimension = (int) (Math.floor(Math.random()*2 )) ; // assume dimension 2
			double oldValue = serum.getParameterValue(whichDimension);
			double change = Math.random()*WALK_SIZE- WALK_SIZE/2 ;
			double value = oldValue+ change;
			
			serum.setParameterValue(whichDimension, value);//WAIT.. IF REJECT, DOES IT RESET?
			


			
			//change I
		//second, find a RANDOM "on" breakpoint and multistep it..
			
			//0. Keep a copy of the original state to calculate the backward move
			int originalNode1 = findAnOnNodeRandomly();

			//System.out.println("Original breakpoint is " + originalNode1 + " and the original AGlikelihood is " + clusterLikelihood.getLogLikelihood());
				
			int[] numStepsFromI_selected =determineTreeNeighborhood(originalNode1, 100000);
			
			//1. Select an unoccupied site within some steps away from it.
			 
			 LinkedList<Integer> possibilities1 = new LinkedList<Integer>();
			 for(int i=0; i < numNodes; i++){
				// System.out.println("#steps from I_selected " + numStepsFromI_selected[i]);
				 //make sure no self select
				 boolean isIn1 = numStepsFromI_selected[i] <= maxNodeLevel && numStepsFromI_selected[i] !=0  && i != rootNum;
				 if(isIn1){
					possibilities1.addLast(new Integer(i));
				 }
			 }//end for
			
			 
			 int numPossibilities1 = possibilities1.size();
			 
			 int whichMove = (int) (Math.floor(Math.random()*numPossibilities1));
				 int site_add1 = possibilities1.get(whichMove).intValue();

			// System.out.println("selected node = " + site_add1 + " that's " + numStepsFromI_selected[site_add1] + " steps from " + originalNode1);
				 
			indicators.setParameterValue(site_add1,  1); //set the new selected index to the new node.
			indicators.setParameterValue(originalNode1, 0); //set the new selected index to the new node.
					

			double[] oldValues = mu.getParameter(site_add1).getParameterValues();
			//double[] oldValues = mu.getParameter(site_add1+1).getParameterValues();	
			//System.out.println(oldValues[0]  + ", " + oldValues[1]); 
			double[] mean = new double[2];
			mean[0] = 0;
			mean[1] = 0;
			double[][] precisionM = new double[2][2];
			//double precision = 1/TreeClusterViruses.getSigmaSq();
			double precision = muPrecision.getParameterValue(0);

			precisionM[0][0] = precision;
			precisionM[0][1] = 0;
			precisionM[1][0] = 0;
			precisionM[1][1] = precision;
			
			
			double[] values = MultivariateNormalDistribution.nextMultivariateNormalPrecision(mean, precisionM);
			//System.out.println(values[0]  + ", " + values[1]); 
			
			mu.getParameter(site_add1).setParameterValue(0,values[0]);
			mu.getParameter(site_add1).setParameterValue(1,values[1]);
			//mu.getParameter(site_add1+1).setParameterValue(0,values[0]);
			//mu.getParameter(site_add1+1).setParameterValue(1,values[1]);

			
			double logHastingRatio = MultivariateNormalDistribution.logPdf(oldValues, mean, precision, 1) - MultivariateNormalDistribution.logPdf(values, mean, precision, 1) ; 
			return logHastingRatio;
	}



	private double Propose_YandI() {
		//first, find a random Y and walk
			int serum_selected = (int) (Math.floor(Math.random()*getNumSera()));
			
			MatrixParameter serumLocations = getSerumLocationsParameter();
			Parameter serum = serumLocations.getParameter(serum_selected);
			int whichDimension = (int) (Math.floor(Math.random()*2 )) ; // assume dimension 2
			double oldValue = serum.getParameterValue(whichDimension);
			double change = Math.random()*WALK_SIZE- WALK_SIZE/2 ;
			double value = oldValue+ change;
			
			serum.setParameterValue(whichDimension, value);//WAIT.. IF REJECT, DOES IT RESET?
			


			int rootNum = treeModel.getRoot().getNumber();
			
		//second, find a RANDOM "on" breakpoint and multistep it..
			int originalNode1 = findAnOnNodeRandomly();			//find an on-node	
			
			//0. Keep a copy of the original state to calculate the backward move


			//System.out.println("Original breakpoint is " + originalNode1 + " and the original AGlikelihood is " + clusterLikelihood.getLogLikelihood());
				
			int[] numStepsFromI_selected =determineTreeNeighborhood(originalNode1, 100000);
			
			//1. Select an unoccupied site within some steps away from it.
			 
			 LinkedList<Integer> possibilities1 = new LinkedList<Integer>();
			 for(int i=0; i < numNodes; i++){
				// System.out.println("#steps from I_selected " + numStepsFromI_selected[i]);
				 //make sure no self select
				 boolean isIn1 = numStepsFromI_selected[i] <= maxNodeLevel && numStepsFromI_selected[i] !=0  && i != rootNum;
				 if(isIn1){
					possibilities1.addLast(new Integer(i));
				 }
			 }//end for
			
			 
			 int numPossibilities1 = possibilities1.size();
			 
			 int whichMove = (int) (Math.floor(Math.random()*numPossibilities1));
				 int site_add1 = possibilities1.get(whichMove).intValue();

			// System.out.println("selected node = " + site_add1 + " that's " + numStepsFromI_selected[site_add1] + " steps from " + originalNode1);
			indicators.setParameterValue(originalNode1, 0); //set the old selected index off
			indicators.setParameterValue(site_add1, 1); //set the new selected index on
					
		
			
						
			
			//it may be more efficient to find a random breakpoint that's closest to it... but it would be hard to code now..
			
			
		
		//what's more important? walk E or mu?
		
		//and move it toegether.
		
		//FOR NOW, don't care about the backward move.. (and MH ratio)..
		//just want to see even if MH ratio is 1, does it ever get accepted..
			return 0;
	}



	private double Propose_YandMu() {
		//first, find a random Y and walk
			int serum_selected = (int) (Math.floor(Math.random()*getNumSera()));
			
			MatrixParameter serumLocations = getSerumLocationsParameter();
			Parameter serum = serumLocations.getParameter(serum_selected);
			int whichDimension = (int) (Math.floor(Math.random()*2 )) ; // assume dimension 2
			double oldValue = serum.getParameterValue(whichDimension);
			double change = Math.random()*WALK_SIZE- WALK_SIZE/2 ;
			double value = oldValue+ change;
			
			serum.setParameterValue(whichDimension, value);//WAIT.. IF REJECT, DOES IT RESET?
			


			
			int selectedIndex = findAnOnNodeRandomly();			//find an on-node	
			
			
			double[] oldValues = mu.getParameter(selectedIndex).getParameterValues();	

			//double[] oldValues = mu.getParameter(selectedIndex+1).getParameterValues();	
			//System.out.println(oldValues[0]  + ", " + oldValues[1]); 
			double[] mean = new double[2];
			mean[0] = 0;
			mean[1] = 0;
			double[][] precisionM = new double[2][2];
			//double precision = 1/TreeClusterViruses.getSigmaSq();
			double precision = muPrecision.getParameterValue(0);

			precisionM[0][0] = precision;
			precisionM[0][1] = 0;
			precisionM[1][0] = 0;
			precisionM[1][1] = precision;
			
			
			double[] values = MultivariateNormalDistribution.nextMultivariateNormalPrecision(mean, precisionM);
			//System.out.println(values[0]  + ", " + values[1]);
			mu.getParameter(selectedIndex).setParameterValue(0,values[0]);
			mu.getParameter(selectedIndex).setParameterValue(1,values[1]);
			//mu.getParameter(selectedIndex+1).setParameterValue(0,values[0]);
			//mu.getParameter(selectedIndex+1).setParameterValue(1,values[1]);

			
			double logHastingRatio = MultivariateNormalDistribution.logPdf(oldValues, mean, precision, 1) - MultivariateNormalDistribution.logPdf(values, mean, precision, 1) ; 
			
			//System.out.println("logHastingRatio = " + logHastingRatio);
			
			
			// but hey, this is not moving the first node.. (it's okay for now)
			
			
			//System.out.println("The first node selected is " + site_add);
			return logHastingRatio;
	}



	private double Proposal_multistepOnNode() {
		
		
		int rootNodeNum = treeModel.getRoot().getNumber();
		//unlike the old version, self-move isn't allowed.
		
		int originalNode1 = findAnOnNodeRandomly();			//find an on-node
		
		//System.out.print("Try " + originalNode1);
		int[] numStepsFromI_selected =determineTreeNeighborhood(originalNode1, 100000);
		
		//1. Select an unoccupied site within some steps away from it.	 
		 LinkedList<Integer> possibilities1 = new LinkedList<Integer>();
		 for(int i=0; i < numNodes; i++){
			// System.out.println("#steps from I_selected " + numStepsFromI_selected[i]);
			 //make sure no self select
			 boolean isIn1 = numStepsFromI_selected[i] <= maxNodeLevel && numStepsFromI_selected[i] !=0 && i != rootNodeNum;
			 if(isIn1){
				possibilities1.addLast(new Integer(i));
			 }
		 }//end for		
		 
		 int numPossibilities1 = possibilities1.size();		 
		 int whichMove = (int) (Math.floor(Math.random()*numPossibilities1)); //choose from possibilities
			 int site_add1 = possibilities1.get(whichMove).intValue();
			 
			curNode = site_add1;

		//	 System.out.println(" and select " + site_add1);
		// System.out.println("selected node = " + site_add1 + " that's " + numStepsFromI_selected[site_add1] + " steps from " + originalNode1);
			 
		indicators.setParameterValue(originalNode1, 0); // the existing indicator is now set to be off
		indicators.setParameterValue(site_add1, 1); //set the new selected index to the new node.
				
	
		//backward calculation
		int[] numStepsBackward =determineTreeNeighborhood(site_add1, 100000);
		
		//1. Select an unoccupied site within some steps away from it.
		 LinkedList<Integer> possibilities2 = new LinkedList<Integer>();
		 for(int i=0; i < numNodes; i++){
			// System.out.println("#steps from I_selected " + numStepsFromI_selected[i]);
			 //make sure no self select
			 boolean isIn2 = numStepsBackward[i] <= maxNodeLevel && numStepsBackward[i] !=0 && i != rootNodeNum;
			 if(isIn2){
				possibilities2.addLast(new Integer(i));
			 }
		 }//end for
		 int numPossibilities2 = possibilities2.size();		 
		 //	 System.out.println("numPossibilities1=" + numPossibilities1 + " numPossibilities2 = " + numPossibilities2);

  		double logHastingRatio = Math.log( (1/ (double)numPossibilities2) / (1/ (double)numPossibilities1)  );
  		//System.out.println("logHastingRatio = " + logHastingRatio);
		return logHastingRatio;
	}

	private double Proposal_flipI(){
		int node = findNodeRandomly();
		if((int)indicators.getParameterValue(node) == 0 ){
			indicators.setParameterValue(node, 1);
		}
		else{
			indicators.setParameterValue(node, 0);
		}
		return(0);
	}

	private double Proposal_flipIandChangeMu() {
		int node = findNodeRandomly();
		if((int)indicators.getParameterValue(node) == 0 ){
			indicators.setParameterValue(node, 1);
		}
		else{
			indicators.setParameterValue(node, 0);
		}
		
		double[] oldValues = mu.getParameter(node).getParameterValues();
		//double[] oldValues = mu.getParameter(node+1).getParameterValues();	
		//System.out.println(oldValues[0]  + ", " + oldValues[1]); 
		double[] mean = new double[2];
		mean[0] = 0;
		mean[1] = 0;
		double[][] precisionM = new double[2][2];
		//double precision = 1/TreeClusterViruses.getSigmaSq();
		double precision = muPrecision.getParameterValue(0);

		precisionM[0][0] = precision;
		precisionM[0][1] = 0;
		precisionM[1][0] = 0;
		precisionM[1][1] = precision;
		
		double[] values = MultivariateNormalDistribution.nextMultivariateNormalPrecision(mean, precisionM);
		//System.out.println(values[0]  + ", " + values[1]); 
		mu.getParameter(node).setParameterValue(0,values[0]);
		mu.getParameter(node).setParameterValue(1,values[1]);
		//mu.getParameter(node+1).setParameterValue(0,values[0]);
		//mu.getParameter(node+1).setParameterValue(1,values[1]);
		
		double logHastingRatio = MultivariateNormalDistribution.logPdf(oldValues, mean, precision, 1) - MultivariateNormalDistribution.logPdf(values, mean, precision, 1) ; 		
		return(logHastingRatio);	
		
		
		
		//int dimSelect = (int) Math.floor( Math.random()* 2 );   		  	    		
		//double change = Math.random()*WALK_SIZE- WALK_SIZE/2 ; 	
		//double originalValue = mu.getParameter(node +1).getParameterValue(dimSelect);				
		//mu.getParameter(node + 1).setParameterValue(dimSelect, originalValue + change);
	}


	private double Proposal_changeRootMuWalk(){
		int dimSelect = (int) Math.floor( Math.random()* 2 );   		  	    		
		double change = Math.random()*WALK_SIZE- WALK_SIZE/2 ;
		
		int rootNum = treeModel.getRoot().getNumber();
	
		double originalValue = mu.getParameter(rootNum).getParameterValue(dimSelect);		
		mu.getParameter(rootNum).setParameterValue(dimSelect, originalValue + change);
		//double originalValue = mu.getParameter(0).getParameterValue(dimSelect);		
		//mu.getParameter(0).setParameterValue(dimSelect, originalValue + change);
		return(0);
	
	}
	
	//Instead of sampling from the prior, I will perform a walk to fine tune things
	private double Proposal_changeAnOnMuWalk() {
		int on_mu =  findAnOnNodeIncludingRootRandomly();

	
		int dimSelect = (int) Math.floor( Math.random()* 2 );   		  	    		
		double change = Math.random()*WALK_SIZE- WALK_SIZE/2 ; 
		double originalValue = mu.getParameter(on_mu).getParameterValue(dimSelect);		
		mu.getParameter(on_mu).setParameterValue(dimSelect, originalValue + change);
		//double originalValue = mu.getParameter(on_mu +1).getParameterValue(dimSelect);		
		//mu.getParameter(on_mu + 1).setParameterValue(dimSelect, originalValue + change);
		return(0);
	}
	
	
	private double Proposal_changeAnOnMuFromPrior(){
		
		int on_mu =  findAnOnNodeIncludingRootRandomly();

	
		double[] oldValues = mu.getParameter(on_mu ).getParameterValues(); //
		
		//double[] oldValues = mu.getParameter(on_mu + 1).getParameterValues();	// this is not +1 because the root's mu can also be changed
		//System.out.println(oldValues[0]  + ", " + oldValues[1]); 
		double[] mean = new double[2];
		//mean[0] = 0;
		mean[0] = muMean.getParameterValue(0);
		mean[1] = 0;
		double[][] precisionM = new double[2][2];
		//double precision = 1/TreeClusterViruses.getSigmaSq();
		double precision = muPrecision.getParameterValue(0);

		precisionM[0][0] = precision;
		precisionM[0][1] = 0;
		precisionM[1][0] = 0;
		precisionM[1][1] = precision;
		
		double[] values = MultivariateNormalDistribution.nextMultivariateNormalPrecision(mean, precisionM);
		mu.getParameter(on_mu ).setParameterValue(0,values[0]);  
		mu.getParameter(on_mu ).setParameterValue(1,values[1]);  
		
		//System.out.println(values[0]  + ", " + values[1]); 
		//mu.getParameter(on_mu + 1).setParameterValue(0,values[0]);  // this is not +1 because the root's mu can also be changed
		//mu.getParameter(on_mu + 1).setParameterValue(1,values[1]);  // this is not +1 because the root's mu can also be changed
		
		double logHastingRatio = MultivariateNormalDistribution.logPdf(oldValues, mean, precision, 1) - MultivariateNormalDistribution.logPdf(values, mean, precision, 1) ; 		
		return(logHastingRatio);
		
	}


	private double Proposal_changeMuWalk() {

		//int nodeSelect = (int) Math.floor( Math.random()* (numNodes + 1) ); 		//pick from index... 0,  to numNodes+1
		int nodeSelect = (int) Math.floor( Math.random()* (numNodes  ) ); 		//pick from index... 0,  to numNodes+1
		
		
		//unbounded walk
		int dimSelect = (int) Math.floor( Math.random()* 2 );   		  	    		
		double change = Math.random()*WALK_SIZE- WALK_SIZE/2 ; 	
		double originalValue = mu.getParameter(nodeSelect).getParameterValue(dimSelect);		
		mu.getParameter(nodeSelect).setParameterValue(dimSelect, originalValue + change);
		return(0);
  				
	}


	private double Proposal_changeMuFromPrior() {

		//int groupSelect = (int) Math.floor( Math.random()* (numNodes + 1) ); 		//pick from index... 0,  to numNodes+1
		int groupSelect = (int) Math.floor( Math.random()* (numNodes ) ); 		//pick from index... 0,  to numNodes+1
		
		double[] oldValues = mu.getParameter(groupSelect).getParameterValues();	// this is not +1 because the root's mu can also be changed
		//System.out.println(oldValues[0]  + ", " + oldValues[1]); 
		double[] mean = new double[2];
		//mean[0] = 0;
		mean[0] = muMean.getParameterValue(0);
		//System.out.println("mean[0] = " + muMean.getParameterValue(0));
		mean[1] = 0;
		double[][] precisionM = new double[2][2];
		//double precision = 1/TreeClusterViruses.getSigmaSq();
		double precision = muPrecision.getParameterValue(0);

		precisionM[0][0] = precision;
		precisionM[0][1] = 0;
		precisionM[1][0] = 0;
		precisionM[1][1] = precision;
		
		double[] values = MultivariateNormalDistribution.nextMultivariateNormalPrecision(mean, precisionM);
		//System.out.println(values[0]  + ", " + values[1]); 
		mu.getParameter(groupSelect).setParameterValue(0,values[0]);  // this is not +1 because the root's mu can also be changed
		mu.getParameter(groupSelect).setParameterValue(1,values[1]);  // this is not +1 because the root's mu can also be changed
		
		double logHastingRatio = MultivariateNormalDistribution.logPdf(oldValues, mean, precision, 1) - MultivariateNormalDistribution.logPdf(values, mean, precision, 1) ;
		
		
		
		//Need this to not screw up the acceptance probability when the indicator is off...
		//if I am using mu_i = 0 | I_i = 0,  the P(mu_i = 1 | I_i = 0) and so this mu_i won't be contributing in the TreeClusterViruses's as a normal prior likelihood
		//and if I am keeping this proposal to latently move the "mu prime", then I don't want selection done on this proposal if
		// I_i = 0.
		//while I think it won't screw up the calculation, it will mess up mixing.
		
		//if((int)indicators.getParameterValue(groupSelect) == 0){
			//logHastingRatio = 0;
		//}
		
		
		return(logHastingRatio);
		
  				
	}



	private double Proposal_changeToAnotherNodeOn() {

			//change another node to be "on" 		
		//System.out.println("a");
	  		//find an on node, turn it off
			int onNode = findAnOnNodeRandomly();
	  		indicators.setParameterValue( onNode , 0);
	  		
	  		//System.out.println("b");
	  		//find an off node, turn it on
	  		int offNode = findAnOffNodeRandomly();
	  		indicators.setParameterValue(offNode  ,1);
	  		
	  		//System.out.println("Change the indicator only - indicator " + I_selected + " to " + site_add);
	  		//System.out.println("Excision point = " + excisionPoints.getParameterValue(I_selected));

	  		//System.out.println("done");
	  		
	  		curNode = offNode;
		
		return(0);
	}


	
	
	

    
    private void test1(){
 		//test whether 
 		//Propose_changeMuAndBalance() and Proposal_changeAnOnMuWalk() are indeed different.
 		
 		//first load initial serum location, mu, and status.
 		System.out.println("Test whether Propose_changeMuAndBalance() and Proposal_changeAnOnMuWalk() are implemented correctly");
 		
 		
 		
	   	System.out.print("  [");
	   	for(int i=0; i < numNodes; i++){
	   		if( (int)indicators.getParameterValue(i) == 1){
	   			System.out.print(i + " ");
	   		}
	   	}
	   	System.out.println("]");
 		
 		Propose_changeMuAndBalance();
 		
 		
 		//Proposal_changeAnOnMuWalk();
 		
 		
 		
 		
 		
 		
 		
 		
 		
 		
 		
 		
 		
 		
 		
 		
 		
 		
 		System.exit(0);
 		
 		
 		//old test
 		
 		
 		int originalNode1 = 605;
 		
 		double originalLikelihood = clusterLikelihood.getLogLikelihood();
 		
 		int[] distance = determineTreeNeighborhood(originalNode1, 5);
 	

		   	
 		for(int x=0; x< 1000; x++){
 			
 			int newNode = 604;

 		indicators.setParameterValue( originalNode1, 0);
 		indicators.setParameterValue( newNode, 1);
 		

 		
 		 		
 		
		//Flip mu - so the neighbor that replaces the original node now also inherits the existing node's mu
 		Parameter originalNodeMu = mu.getParameter(originalNode1); //offset of 1
		//Parameter originalNodeMu = mu.getParameter(originalNode1+1); //offset of 1
		double[] tmp = originalNodeMu.getParameterValues();

		Parameter newMu = mu.getParameter(newNode); //offset of 1
		//Parameter newMu = mu.getParameter(newNode+1); //offset of 1
		double[] tmpNew = newMu.getParameterValues();
		
		double change0 =  Math.random()*4 - 2;
		double change1 =  Math.random()*4 - 2;
		
		originalNodeMu.setParameterValue(0, tmpNew[0]);
		originalNodeMu.setParameterValue(1, tmpNew[1]);

		newMu.setParameterValue(0, tmp[0] + change0);
		newMu.setParameterValue(1, tmp[1] + change1);	 		

 		 
				
		//1. Update the cluster labels, after the breakpoints and status parameters may have changed.
	//		setClusterLabelsUsingIndicators();
			//setClusterLabelsArray(newClusterLabelArray);	
			//relabelClusterLabelsArray(newClusterLabelArray, oldClusterLabelArray);
			//convertClusterLabelsArrayToParameter( newClusterLabelArray);
			//oldClusterLabelArray = newClusterLabelArray; //the oldClusterLabelArray gets the current labels, so next time this is updated.
		

		//2. Update the virus locations (and offsets), given ...
		TreeClusteringSharedRoutines.updateUndriftedVirusLocations(numNodes, numdata, treeModel, virusLocationsTreeNode, indicators, mu, virusLocations, correspondingTreeIndexForVirus);
			//setVirusLocationAutoCorrelatedModel(); //set virus locations, given the breakpoints,status, and mu parameters
		
			double testLikelihood = clusterLikelihood.getLogLikelihood();
			double diff =  testLikelihood - originalLikelihood ;
			
			if(diff > 0){
				System.out.print("***");
			}
		System.out.println("logL= " + testLikelihood + " and diff = " + diff);
 		
 		
	
	   	
	   	
	   	

		//revert back
 		indicators.setParameterValue( originalNode1, 1);
 		indicators.setParameterValue( newNode, 0);

		originalNodeMu.setParameterValue(0, tmp[0]);
		originalNodeMu.setParameterValue(1, tmp[1]);

		newMu.setParameterValue(0, tmpNew[0]);
		newMu.setParameterValue(1, tmpNew[1]);	 		


	   	
 			
 		}
	   	
	   	System.exit(0);


 		
    }
    
    private void test2(){

 	   	System.out.print("  [");
	   	for(int i=0; i < numNodes; i++){
	   		if( (int)indicators.getParameterValue(i) == 1){
	   			System.out.print(i + " ");
	   		}
	   	}
	   	System.out.println("]");
		
		
	   	
	  	//indicators.setParameterValue( 436, 0);
	   	//indicators.setParameterValue( 549, 0);
	//   	indicators.setParameterValue( 615, 0);   
	   	//indicators.setParameterValue(648,0);
	   	//indicators.setParameterValue(673,0);
	   	//indicators.setParameterValue(794,0);
	   	//indicators.setParameterValue(785,0);
	   	//indicators.setParameterValue(690,0);
	   	
	   	//Note: if 615 is not turned on, 604 would be much superior to 605.
	   	
 	//since I changed indicators.. should do this before calculating the original likelihood
 		

		//1. Update the cluster labels, after the breakpoints and status parameters may have changed.
	 //  		setClusterLabelsUsingIndicators();
	   		//setClusterLabelsArray(newClusterLabelArray);	
			//relabelClusterLabelsArray(newClusterLabelArray, oldClusterLabelArray);
			//convertClusterLabelsArrayToParameter( newClusterLabelArray);
			//oldClusterLabelArray = newClusterLabelArray; //the oldClusterLabelArray gets the current labels, so next time this is updated.
		

		//2. Update the virus locations (and offsets), given ...
	   	TreeClusteringSharedRoutines.updateUndriftedVirusLocations(numNodes, numdata, treeModel, virusLocationsTreeNode, indicators, mu, virusLocations, correspondingTreeIndexForVirus);
			//setVirusLocationAutoCorrelatedModel(); //set virus locations, given the breakpoints,status, and mu parameters
		
			
			
			double originalLikelihood = clusterLikelihood.getLogLikelihood();
			
			System.out.println("originalLikelihood = " + originalLikelihood);
 		
 		
 		
 		int originalNode1 = 605;
 		
 		
 		
 		int[] distance = determineTreeNeighborhood(originalNode1, 5);
 		
 		
 		for(int g=0; g < distance.length; g++){
 			if(distance[g] < 5 && distance[g] >0){
 		 		int newNode = g;
 				System.out.print(g + " distance=" + distance[g] + "\t");
 				
 				
	 		if( (int) indicators.getParameterValue(newNode) == 1){
	 			System.out.print("Node already on!!!\n");
	 		}
	 		else{

 		indicators.setParameterValue( originalNode1, 0);
 		indicators.setParameterValue( newNode, 1);
 		

 		
 		 		
 		
		//Flip mu - so the neighbor that replaces the original node now also inherits the existing node's mu
		//Parameter originalNodeMu = mu.getParameter(originalNode1+1); //offset of 1
 		Parameter originalNodeMu = mu.getParameter(originalNode1); 
		double[] tmp = originalNodeMu.getParameterValues();

		//Parameter newMu = mu.getParameter(newNode+1); //offset of 1
		Parameter newMu = mu.getParameter(newNode); 
		double[] tmpNew = newMu.getParameterValues();
		
		originalNodeMu.setParameterValue(0, tmpNew[0]);
		originalNodeMu.setParameterValue(1, tmpNew[1]);

		newMu.setParameterValue(0, tmp[0]);
		newMu.setParameterValue(1, tmp[1]);	 		

 		 
				
		//1. Update the cluster labels, after the breakpoints and status parameters may have changed.
	//		setClusterLabelsUsingIndicators();
			//setClusterLabelsArray(newClusterLabelArray);	
			//relabelClusterLabelsArray(newClusterLabelArray, oldClusterLabelArray);
			//convertClusterLabelsArrayToParameter( newClusterLabelArray);
			//oldClusterLabelArray = newClusterLabelArray; //the oldClusterLabelArray gets the current labels, so next time this is updated.
		

		//2. Update the virus locations (and offsets), given ...
		TreeClusteringSharedRoutines.updateUndriftedVirusLocations(numNodes, numdata, treeModel, virusLocationsTreeNode, indicators, mu, virusLocations, correspondingTreeIndexForVirus);
			//setVirusLocationAutoCorrelatedModel(); //set virus locations, given the breakpoints,status, and mu parameters
		
			double testLikelihood = clusterLikelihood.getLogLikelihood();
			double diff =  testLikelihood - originalLikelihood ;
			
			if(diff > -10){
				System.out.print("***");
			}
		System.out.println("logL= " + testLikelihood + " and diff = " + diff);
 		
 		
	  // 	System.out.print("  [");
	   //	for(int i=0; i < numNodes; i++){
	   //		if( (int)indicators.getParameterValue(i) == 1){
	   //			System.out.print(i + " ");
	   //		}
	   	//}
	   	//System.out.println("]");
	   	
	   	
	   	

		//revert back
 		indicators.setParameterValue( originalNode1, 1);
 		indicators.setParameterValue( newNode, 0);

		originalNodeMu.setParameterValue(0, tmp[0]);
		originalNodeMu.setParameterValue(1, tmp[1]);

		newMu.setParameterValue(0, tmpNew[0]);
		newMu.setParameterValue(1, tmpNew[1]);	 		
	 		}

	   	
 			}
 		}
	   	
	   	System.exit(0);

    }
    
    
    private void test3(){
 		System.out.println("Turn a new node on");

 		//696 and 0.1 ,0 works
	   	//manually flip a new node on
 		for(int newNode = 0; newNode < 803; newNode++){
 			
				System.out.print(newNode +  "\t");

 		if( (int) indicators.getParameterValue(newNode) == 1){
 			System.out.print("Node already on!!!\n");
 		}
 		else{
 			
 			
	 		double originalLikelihood = clusterLikelihood.getLogLikelihood();

 			indicators.setParameterValue( newNode, 1);
			//Parameter newMu = mu.getParameter(newNode+1); //offset of 1
 			Parameter newMu = mu.getParameter(newNode);
			double[] tmpNew = newMu.getParameterValues();
			newMu.setParameterValue(0, 0);
			newMu.setParameterValue(1, 0);
		
	 							
			//1. Update the cluster labels, after the breakpoints and status parameters may have changed.
		//		setClusterLabelsUsingIndicators();
				//setClusterLabelsArray(newClusterLabelArray);	
				//relabelClusterLabelsArray(newClusterLabelArray, oldClusterLabelArray);
				//convertClusterLabelsArrayToParameter( newClusterLabelArray);
				//oldClusterLabelArray = newClusterLabelArray; //the oldClusterLabelArray gets the current labels, so next time this is updated.
			

			//2. Update the virus locations (and offsets), given ...
			TreeClusteringSharedRoutines.updateUndriftedVirusLocations(numNodes, numdata, treeModel, virusLocationsTreeNode, indicators, mu, virusLocations, correspondingTreeIndexForVirus);
				//setVirusLocationAutoCorrelatedModel(); //set virus locations, given the breakpoints,status, and mu parameters
			
				double testLikelihood = clusterLikelihood.getLogLikelihood();
				double diff =  testLikelihood - originalLikelihood ;
				
				if(diff > 0){
					System.out.print("***");
				}
			System.out.println("logL= " + testLikelihood + " and diff = " + diff);
	 		
			  // 	System.out.print("  [");
			   //	for(int i=0; i < numNodes; i++){
			   //		if( (int)indicators.getParameterValue(i) == 1){
			   //			System.out.print(i + " ");
			   //		}
			   	//}
			   //	System.out.println("]");
			   	
			   	//revert back
		 		indicators.setParameterValue( newNode, 0);
				newMu.setParameterValue(0, tmpNew[0]);
				newMu.setParameterValue(1, tmpNew[1]);	
			
 		}
 		
 		} //for
 		
 		System.exit(0);
    }

	
	
	
	
	
	
	
	
	//===============================================================================================
	//===============================================================================================
	
	//  BELOW IS A LIST OF HELPER ROUTINES
	
	//===============================================================================================
	//===============================================================================================
	
	
	public double getNumSera() {
		return serumLocations.getParameterCount();
	}

	public MatrixParameter getSerumLocationsParameter() {
		return serumLocations;
	}
  
	
	
	private LinkedList<Integer> findActiveBreakpointsChildren(int selectedNodeNumber) {
		
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



	private int checkSiteHasBeenAddedToOnIndicators(int curTest){
		
		int hasBeenAdded=0;
		if((int) indicators.getParameterValue(curTest) == 1){
			hasBeenAdded=1;
		}
		
		return(  hasBeenAdded );
	}


	//for Gibbs move
    private double[] calculateConditionalDistribution(int index) {
		double []logNumeratorProb = new double[numNodes];

		//calculate the distribution for calculating introducing an excision point in each node
		for(int curTest=0; curTest < numNodes; curTest++){	
			
			int hasBeenAdded = checkSiteHasBeenAddedToOnIndicators(curTest); 			//check if a site has already been added
  			if(hasBeenAdded ==0){

  		    	indicators.setParameterValue(curTest,1);
  		    	
  		    	updateClusterLabelsAndVirusLocationsGivenBreakPointsAndStatus();	    		  					  				
   				logNumeratorProb[curTest] = clusterLikelihood.getLogLikelihood(); 	//Calculate likelihood
   				
   				
		    	indicators.setParameterValue(curTest,0); //set back to original
			    
  			}
  			else{
	  			logNumeratorProb[curTest]  = Double.NEGATIVE_INFINITY; //dummy probability
  			}
  			

		} //finished curTest
		
		 double []condDistribution = calculateConditionalProbabilityGivenLogNumeratorProb(logNumeratorProb);
		
//		 System.out.println("-----");
//		 for(int i=0; i < numNodes; i++){
//			 if(condDistribution[i] > 0.0000001){
//				 System.out.println("node " + i + " p=" + condDistribution[i]);
//			 }
//		 }
//		 System.out.println("-----");
		 
		 
		 updateClusterLabelsAndVirusLocationsGivenBreakPointsAndStatus();	
		 
		 return condDistribution;
	}

    
    
	private void updateClusterLabelsAndVirusLocationsGivenBreakPointsAndStatus() {

		
			updateClusterLabelsWhileKeepingLablesConsistent();	  //CURRENTLY DOES NOT WORK..				
			
			//setVirusLocationAndOffsets();  //this uses the clusterLabels parameter
			TreeClusteringSharedRoutines.updateUndriftedVirusLocations(numNodes, numdata, treeModel, virusLocationsTreeNode, indicators, mu, virusLocations, correspondingTreeIndexForVirus);
			//setVirusLocationAutoCorrelatedModel(); //which depends on the status and breakpoints
					
	}



	private void updateClusterLabelsWhileKeepingLablesConsistent() {

		/*
		int old
    	//use the tree to re-partition according to the change.
		clusterLabelArray = setClusterLabelsByTestCutNodeByNodeOrder(testCutNode); //note that instead of using the indicators, it uses the testCutNodes directly
		relabelClusterLabels(clusterLabelArray, oldclusterLabelArray); //will move it out
		
		//set cluster label parameter for testing 					
		for(int i=0; i < numdata; i++){
			clusterLabels.setParameterValue(i, clusterLabelArray[i]);
		}
		*/
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
	
	


	//RETIRE
	/* no longer needed, I think
	private int findAnUnoccupiedSite() {
    	
    	int hasBeenAdded = 1;
    	int site_add = -1;
    	while(hasBeenAdded==1){
  			site_add = (int) Math.floor( Math.random()*numNodes );
  			hasBeenAdded=0;
 			if( (int) indicators.getParameterValue(site_add) == 1){
 				hasBeenAdded=1;
 				break;
 			}
    	}
	  	
		return site_add;
	}
	*/


	
	//may be very inefficient
	private int findNodeRandomly() {
    	int nodeIndex= 0;
    	int I_selected = -1;
  		while(nodeIndex ==0){
  			I_selected = (int) (Math.floor(Math.random()*numNodes));
  			
  			if(I_selected == treeModel.getRoot().getNumber()){
  				nodeIndex = 0;
  			}
  			else{
  	  			nodeIndex = 1;  				
  			}
  		}    	  		
    	return I_selected;
	}

	
	//may be very inefficient
	private int findAnOnNodeIncludingRootRandomly() {
    	int isOn= 0;
    	int I_selected = -1;
  		while(isOn ==0){
  			I_selected = (int) (Math.floor(Math.random()*numNodes));
  			isOn = (int) indicators.getParameterValue(I_selected);  			
  		}    	  		
  		
    	return I_selected;
	}


	//may be very inefficient
	private int findAnOnNodeRandomly() {
    	int isOn= 0;
    	int I_selected = -1;
  		while(isOn ==0){
  			I_selected = (int) (Math.floor(Math.random()*numNodes));
  			isOn = (int) indicators.getParameterValue(I_selected);
  			
  			if(I_selected == treeModel.getRoot().getNumber()){
  				isOn = 0;
  			}
  			
  		}    	  		
  		
    	return I_selected;
	}

	
	private int findAnOffNodeRandomly() {
    	int isOn= 1;
    	int I_selected = -1;
  		while(isOn ==1){
  			I_selected = (int) (Math.floor(Math.random()*numNodes));
  			isOn = (int) indicators.getParameterValue(I_selected);
  		}    	  		
  		
    	return I_selected;
	}

	
	
	

	/*
	private void updateK() {

    	//K is changed accordingly..
		int K_count = 0; //K_int gets updated
		for(int i=0; i < numNodes; i++){
			K_count += (int) indicators.getParameterValue(i);
		}
		//System.out.println("K now becomes " + K_count);
		K.setParameterValue(0, K_count); //update   
 							
	}
*/


/*
	private void convertClusterLabelsArrayToParameter(int[] clusterLabel){
    	for(int i=0; i < clusterLabel.length; i++){
    		clusterLabels.setParameterValue(i, clusterLabel[i]);
    	}
    }
  */
	
	private void relabelClusterLabelsArray(int[] clusterLabel, int[] oldclusterLabel) {

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
    				
    				if( clusterLabel[i] != oldclusterLabel[i]){
    					System.out.println("conversion occurred");
    				}
    			}
    			else{
    				maxOldLabel++;
    				m.put(new Integer(clusterLabel[i]), new Integer(maxOldLabel));
    			}
    			
    		}

    		clusterLabel[i] = m.get(new Integer( clusterLabel[i])).intValue();
    		
    	}
	}

    
    
    
/*
	private void setMembershipTreeToVirusIndexes(){

  	   //I suspect this is an expensive operation, so I don't want to do it many times,
  	   //which is also unnecessary  - MAY have to update whenever a different tree is used.
         correspondingTreeIndexForVirus = new int[numdata]; 
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
    }
*/

	

	private void PrintsetMembershipTreeToVirusIndexes(){

  	   //I suspect this is an expensive operation, so I don't want to do it many times,
  	   //which is also unnecessary  - MAY have to update whenever a different tree is used.
         correspondingTreeIndexForVirus = new int[numdata]; 
         for(int i=0; i < numdata; i++){
  		   Parameter v = virusLocations.getParameter(i);
  		   String curName = v.getParameterName();
  		   System.out.print(curName);
  		   int isFound = 0;
      	   for(int j=0; j < numNodes; j++){
      		   String treeId = treeModel.getTaxonId(j);
      		   if(curName.equals(treeId) ){
      			   System.out.print("  isFound at j=" + j);
      			   correspondingTreeIndexForVirus[i] = j;
      			   System.out.println(" has clusterLabel = " + clusterLabelsTreeNode.getParameterValue(j));
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


	
	
	
	
	
	
/*
	//Obsolete
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
*/

    
/*
	private void setVirusLocationAutoCorrelatedModel() {
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
	*/



	private int[] determineTreeNeighborhood(int curElementNumber, int maxDepth ) {
			int numNodes = treeModel.getNodeCount();
 			
			//Determining the number of steps from the original site
			int []numStepsFromOrigin = new int[numNodes];
			for(int i=0; i < numNodes; i++){
				numStepsFromOrigin[i] = 100000;
			}
 		
  		//System.out.println("Excision point = " + excisionPoints.getParameterValue(I_selected));

  		
  			//int curElementNumber =(int) indicators.getParameterValue(I_selected);
  			int rootElementNumber = curElementNumber;
  			//System.out.println("curElementNumber=" + curElementNumber);
  			NodeRef curElement = treeModel.getNode(curElementNumber); 
  			
  		    LinkedList<NodeRef> visitlist = new LinkedList<NodeRef>();
  		    LinkedList<NodeRef> fromlist = new LinkedList<NodeRef>();
  		    LinkedList<Integer> nodeLevel = new LinkedList<Integer>();
  		    
  		    //LinkedList<Integer> possibilities = new LinkedList<Integer>();
  		    
  		    NodeRef dummyNode = null;
  		    visitlist.add(curElement);
  		    fromlist.add(dummyNode);
  		    nodeLevel.add(new Integer(0));
  		    
  		    //int numVisited = 0;
  		    
  		  //System.out.println("root node " + curElement.getNumber());
		    while(visitlist.size() > 0){
		    	//numVisited++;
		    	
  			if(treeModel.getParent(curElement) != null){
  				//add parent
		  			NodeRef node= treeModel.getParent(curElement);	  		  			
  				if(fromlist.getFirst() != node){
  					if( nodeLevel.getFirst() < maxDepth){
  						visitlist.add(node);
  		  				fromlist.add(curElement);
  		  				nodeLevel.add(new Integer(nodeLevel.getFirst()+1));
  		 // 				System.out.println("node " +  node.getNumber() + " added, parent of " + curElement.getNumber());
  					}
  				}
  			}

			
  			for(int childNum=0; childNum < treeModel.getChildCount(curElement); childNum++){
  				NodeRef node= treeModel.getChild(curElement,childNum);
  				if(fromlist.getFirst() != node){
  					if( nodeLevel.getFirst() < maxDepth){
  						visitlist.add(node);
  						fromlist.add(curElement);
  						nodeLevel.add(new Integer(nodeLevel.getFirst()+1));
  			//			System.out.println("node " +  node.getNumber() + " added, child of " + curElement.getNumber());
  					}
  				}
  	        }
  			
  			
				numStepsFromOrigin[curElement.getNumber()] = nodeLevel.getFirst();

  			
	  			visitlist.pop();
	  			fromlist.pop();
	  			nodeLevel.pop();

  			if(visitlist.size() > 0){
  				curElement = visitlist.getFirst();
  			}
  			
  			
		}
		    
		 			    
		    
		    return(numStepsFromOrigin);
		
	}



	
	private void setClusterLabelsUsingIndicators(){

        int []membership = determine_membership_v2(treeModel);
        
        for(int i=0; i < numdata; i++){    
        	clusterLabels.setParameterValue(i,membership[correspondingTreeIndexForVirus[i]] );
        }
	}
	
	
	
	private void setClusterLabelsTreeNodesUsingIndicators(){
		
        int []membership = determine_membership_v2(treeModel);
        for(int i=0; i < numNodes; i++){
        	clusterLabelsTreeNode.setParameterValue(i, membership[i]);
        }
	}
	
	
	//composite:

	private void CompositeSetClusterLabelsTreeNodesAndVirusesUsingIndicators(){
		//setMembershipTreeToVirusIndexes(); //note: I have to add this in to fix the inconsistency between
		//the clusterLabelsTreeNode and clusterLabels.. 
		//do I really need to do this everytime?
		//I always thought if I use the same tree, it won't change?
		
        int []membership = determine_membership_v2(treeModel);
        for(int i=0; i < numNodes; i++){
        	clusterLabelsTreeNode.setParameterValue(i, membership[i]);
        }
        for(int i=0; i < numdata; i++){    
        	clusterLabels.setParameterValue(i,membership[correspondingTreeIndexForVirus[i]] );
        }
	}

	

    //traverse down the tree, top down, do calculation
     int[] determine_membership_v2(TreeModel treeModel){
	    	
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
	    		if( (int) indicators.getParameterValue(curElement.getNumber() ) == 1) {
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

	
	
	
	
	private void setClusterLabelsArray(int[] clusterLabelArray) {
		int K_int = 0;
        for(int i=0; i < numNodes; i++){
      	   if( (int) indicators.getParameterValue( i ) ==1 ){
      		  K_int++; 
      	   }
        }
       // System.out.println("K_int=" + K_int);
        int numNodes = treeModel.getNodeCount();
        int[] cutNodes = new int[K_int];
 	   int cutNum = 0;
 	   String content = "";
        for(int i=0; i < numNodes; i++){
     	   if( (int) indicators.getParameterValue( i ) == 1 ){
     		   cutNodes[cutNum] = i;
     		   content +=  i +  ",";
     		   cutNum++;
     	   }
     	  
        }
        //System.out.println( content);        


        int []membership = determine_membership(treeModel, cutNodes, K_int);
        
        for(int i=0; i < numdata; i++){    
        	clusterLabelArray[i] = membership[correspondingTreeIndexForVirus[i]];
        }
    	
	}
	
	/*
	private void setClusterLabelsParameter() {
		int K_int = 0;
        for(int i=0; i < numNodes; i++){
      	   if( (int) indicators.getParameterValue( i ) ==1 ){
      		  K_int++; 
      	   }
        }
        int numNodes = treeModel.getNodeCount();
        int[] cutNodes = new int[K_int];
 	   int cutNum = 0;
 	   String content = "";
        for(int i=0; i < numNodes; i++){
     	   if( (int) indicators.getParameterValue( i ) ==1 ){
     		   cutNodes[cutNum] = i;
     		   content +=  i + ",";
     		   cutNum++;
     	   }
     	  
        }
        


        int []membership = determine_membership(treeModel, cutNodes, K_int);
        
        for(int i=0; i < numdata; i++){     	   
     	   clusterLabels.setParameterValue( i, membership[membershipToClusterLabelIndexes[i]]);      	   //The assumption that the first nodes being external node corresponding to the cluster labels IS FALSE, so I have to search for the matching indexes
     	   //Parameter vloc = virusLocations.getParameter(i);
     	   //System.out.println(vloc.getParameterName() + " i="+ i + " membership=" + (int) clusterLabels.getParameterValue(i));
        }
        

    	
	}
*/
    
	    
	
	
    
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
     int[] determine_membership(TreeModel treeModel, int[] cutNodes, int numCuts){
	    	
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
	   

	  
	  public void hotNodeProcedure(){


	    	//update... Actually, I don't think it should be here
	    	if( moveCounter % updateHotNodeFrequencey == 0 ){
	    		System.out.print("Update hot nodes: ");
	    		
	    		for(int i=0; i < numNodes; i++){
	    			if(freqAcceptNode[i] > 0){
	    				hotNodes[i] = 1;
	    				System.out.print(i + " ");
	    			}
	    			else{
	    				hotNodes[i] = 0;
	    			}
	    		}
	    		System.out.println("");
	    		//reset count
	    		for(int i=0; i < numNodes; i++){
	    			freqAcceptNode[i] = 0;
	    		}
	    	}

	  }
	  
	  public void printAcceptance(){
		  	if(moveCounter > BURN_IN  &&  moveCounter % frequencyPrintAcceptance == 0 ){
	        	System.out.println("======================================================");
	        	System.out.println("#\tProposal\tAcceptance Rate");
	        	for(int i=0; i < operatorWeight.length; i++){
	        		if(operatorWeight[i] > 0 ){
	        			System.out.println(i +"\t" + operatorName[i] + "\t" + (double) acceptNum[i] / (double) (acceptNum[i] + (double) rejectNum[i])  + "\taccept=" + acceptNum[i] + " reject=" + rejectNum[i]);
	        		}
	        	}
	        	System.out.println("======================================================");

	        	
	    	  	//reset acceptance
	        	for(int i=0; i < operatorWeight.length; i++){
	        		if(operatorWeight[i] > 0 ){
	        			acceptNum[i]= 0;
	        			rejectNum[i] = 0;
	        		}
	        	}
		  
	        	
	        	/*
	        	System.out.println("Acceptance of flipIBalance by distance:");
	        	for(int i=0; i < 20; i++){
        			System.out.println( ((double)i ) +"\t" + (double) acceptDistance[i] / (double) (acceptDistance[i] + (double) rejectDistance[i])  + "\taccept=" + acceptDistance[i] + " reject=" + rejectDistance[i]);
	        	}
	        	
	        	
	        	for(int i=0; i < 100; i++){
	        		acceptDistance[i] = 0;
	        		rejectDistance[i] = 0;
	        	}
	        	*/
	    	}
		  	
	
	  
	  }
	  
	  /*
	  public void printActiveNodes(){

	    	if( moveCounter % frequencyPrintActive == 0 ){
		   	System.out.print("  [");
		   	for(int i=0; i < numNodes; i++){
		   		if( (int)indicators.getParameterValue(i) == 1){
		   			System.out.print(i + " ");
		   		}
		   	}
		   	System.out.println("]");
	    	}
	  }
	  */
	public void accept(double deviation) {
    	super.accept(deviation);

    	
    	//if(curNode != -1){
    	//	freqAcceptNode[curNode]++;
    	//}
    	//hotNodeProcedure();
    	
    	
    	if(moveCounter > BURN_IN){
    		acceptNum[operatorSelect]++;
    	}
    	printAcceptance();
    	
    	//printActiveNodes();
  
   	   	

        			
    			/*
    			 *     //obsolete - to see which multistep gets acccepted
    		if(operatorSelect ==2){
    			acceptNumStep[ howmanyStepsMultiStep]++;
    			
    			if( howmanyStepsMultiStep > 0){
    	    		System.out.println("accept operator " + operatorSelect+" with overall % accept = " + acceptNum[operatorSelect]/(acceptNum[operatorSelect] + rejectNum[operatorSelect]));

    				System.out.println("   > with step = " + howmanyStepsMultiStep  + " (now only print this if step >0)");
    	  			System.out.print(" dist=[" );
        			for(int i=0; i <= maxNodeLevel; i++){
        				System.out.print(acceptNumStep[i]/(acceptNumStep[i] + rejectNumStep[i]) + ",\t ");
        			}
        			System.out.println("]");
    			}
  
    		}
    		*/
    	
    		/*
          	if(operatorSelect == 15){
          		acceptDistance[ (int)  (muDistance)]++;
          	}
          	*/
    }
    
    public void reject(){
    	super.reject();
    	
    	//hotNodeProcedure();
    	
    	if(moveCounter > BURN_IN){
    		rejectNum[operatorSelect]++;
    	}
    	printAcceptance();
    	
    	//printActiveNodes();
    	
    //obsolete - to see which multistep gets acccepted	
	//	if(operatorSelect ==2){
	//		rejectNumStep[ howmanyStepsMultiStep]++;
	//	}
    	/*
      	if(operatorSelect == 15){
      		rejectDistance[ (int) (muDistance)]++;
      	}
      	*/

    }
    
	


              
            //MCMCOperator INTERFACE
            public final String getOperatorName() {
                return TREE_CLUSTERALGORITHM_OPERATOR;
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
                public final static String SERUMLOCATIONS = "serumLocations";

            	public final static String  MU = "mu";
            	public final static String CLUSTERLABELS = "clusterLabels";
            	//public final static String K = "k";
            	public final static String OFFSETS = "offsets";
      //     	public final static String LOCATION_DRIFT = "locationDrift"; //no longer need
            	
                public final static String CLUSTER_OFFSETS = "clusterOffsetsParameter";
                
            	public final static String INDICATORS = "indicators";

                public final static String EXCISION_POINTS = "excisionPoints";

                public final static String MUPRECISION = "muPrecision";

                
                public final static String FILE_NAME = "fileName";

                public final static String CLUSTERLABELSTREENODE = "clusterLabelsTreeNodes";
                public final static String VIRUSLOCATIONSTREENODE = "virusLocationsTreeNodes";
                
                
                public final static String MU1SCALE = "mu1Scale";
                public final static String MU2SCALE = "mu2Scale";
                public final static String MUMEAN = "muMean";
                
                public String getParserName() {
                    return TREE_CLUSTERALGORITHM_OPERATOR;
                }

                /* (non-Javadoc)
                 * @see dr.xml.AbstractXMLObjectParser#parseXMLObject(dr.xml.XMLObject)
                 */
                public Object parseXMLObject(XMLObject xo) throws XMLParseException {

                	
                    String fileName = xo.getStringAttribute(FILE_NAME);
                    DataTable<String[]> proposalWeightTable;
                    try {
                    	proposalWeightTable = DataTable.Text.parse(new FileReader(fileName), false, false);
                    	
                    	                    	
                    } catch (IOException e) {
                        throw new XMLParseException("Unable to read proposal weight from file: " + e.getMessage());
                    }
                    System.out.println("Loaded proposal weight table file: " + fileName);

                    
                	//System.out.println("Parser run. Exit now");
                	//System.exit(0);

                    double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);

                    
                    XMLObject cxo = xo.getChild(VIRUSLOCATIONS);
                        MatrixParameter virusLocations = (MatrixParameter) cxo.getChild(MatrixParameter.class);
                        
                        
                        cxo = xo.getChild(VIRUSLOCATIONSTREENODE);
                        MatrixParameter virusLocationsTreeNode = (MatrixParameter) cxo.getChild(MatrixParameter.class);

                        cxo = xo.getChild(SERUMLOCATIONS);
                        MatrixParameter serumLocations = (MatrixParameter) cxo.getChild(MatrixParameter.class);

                        
                        cxo = xo.getChild(MU);
                        MatrixParameter mu = (MatrixParameter) cxo.getChild(MatrixParameter.class);

                        cxo = xo.getChild(CLUSTERLABELSTREENODE);
                        Parameter clusterLabelsTreeNode = (Parameter) cxo.getChild(Parameter.class);
                        
                        cxo = xo.getChild(CLUSTERLABELS);
                        Parameter clusterLabels = (Parameter) cxo.getChild(Parameter.class);

                        //cxo = xo.getChild(K); //to be deleted
                        //Parameter k = (Parameter) cxo.getChild(Parameter.class); //to be deleted
                        
                        
          //              cxo = xo.getChild(OFFSETS);
           //             Parameter offsets = (Parameter) cxo.getChild(Parameter.class);
                        Parameter offsets = null;
                        
//                        cxo = xo.getChild(LOCATION_DRIFT);
//                        Parameter locationDrift = (Parameter) cxo.getChild(Parameter.class);
                        
                        Parameter clusterOffsetsParameter = null;
                //        if (xo.hasChildNamed(CLUSTER_OFFSETS)) {
                 //       	clusterOffsetsParameter = (Parameter) xo.getElementFirstChild(CLUSTER_OFFSETS);
                  //      }

                        cxo = xo.getChild(INDICATORS);
                        Parameter indicators = (Parameter) cxo.getChild(Parameter.class);
                        
                        cxo = xo.getChild(MUPRECISION);
                        Parameter muPrecision = (Parameter) cxo.getChild(Parameter.class);
                      
                      
                        
                        cxo = xo.getChild(MU1SCALE);
                        Parameter mu1Scale = (Parameter) cxo.getChild(Parameter.class);
                        
                        cxo = xo.getChild(MU2SCALE);
                        Parameter mu2Scale = (Parameter) cxo.getChild(Parameter.class);
                        
                    TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
                    	
                    	//AGLikelihoodTreeCluster agLikelihood = (AGLikelihoodTreeCluster) xo.getChild(AGLikelihoodTreeCluster.class);
                    AntigenicLikelihood agLikelihood = (AntigenicLikelihood) xo.getChild(AntigenicLikelihood.class);

                    	cxo = xo.getChild(MUMEAN);
                    	Parameter muMean = (Parameter) cxo.getChild(Parameter.class);
                        
                    //return new ClusterAlgorithmOperator(virusLocations, mu, clusterLabels, k, weight, offsets, locationDrift, clusterOffsetsParameter);
                  //      return new TreeClusterAlgorithmOperator(virusLocations, serumLocations, mu, clusterLabels, k, weight, offsets,  clusterOffsetsParameter, indicators, treeModel, agLikelihood);
                    return new TreeClusterAlgorithmOperator(virusLocations, virusLocationsTreeNode, serumLocations, mu,  clusterLabels, weight,  indicators, treeModel, agLikelihood, muPrecision, proposalWeightTable, clusterLabelsTreeNode, mu1Scale, mu2Scale, muMean);
                    

                }

                //************************************************************************
                // AbstractXMLObjectParser implementation
                //************************************************************************

                public String getParserDescription() {
                    return "tree cluster algorithm's main operator.";
                }

                public Class getReturnType() {
                    return TreeClusterAlgorithmOperator.class;
                }


                public XMLSyntaxRule[] getSyntaxRules() {
                    return rules;
                }

                private final XMLSyntaxRule[] rules = {
                        AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
                        
                        new ElementRule(VIRUSLOCATIONS, Parameter.class),
                        new ElementRule(SERUMLOCATIONS, Parameter.class),
                        new ElementRule(MU, Parameter.class),
                        new ElementRule(CLUSTERLABELS, Parameter.class),
                        //new ElementRule(K, Parameter.class), //to be deleted
                  //      new ElementRule(OFFSETS, Parameter.class),
                  //      new ElementRule(LOCATION_DRIFT, Parameter.class), //no longer needed
   //                    
                  //     new ElementRule(CLUSTER_OFFSETS, Parameter.class, "Parameter of cluster offsets of all virus"),  // no longer REQUIRED
                       new ElementRule(INDICATORS, Parameter.class),
                       new ElementRule(EXCISION_POINTS, Parameter.class),
                       new ElementRule(TreeModel.class),
                       new ElementRule(MUPRECISION, Parameter.class),
                       new ElementRule(CLUSTERLABELSTREENODE, Parameter.class),
                       new ElementRule(MU1SCALE, Parameter.class),
                       new ElementRule(MU2SCALE, Parameter.class),
                       new ElementRule(MUMEAN, Parameter.class),
                       new ElementRule(VIRUSLOCATIONSTREENODE, MatrixParameter.class),
                       AttributeRule.newStringRule(FILE_NAME, false, "The name of the file containing the assay table"),
        
            };
            
            };


        
            public int getStepCount() {
                return 1;
            }

        }



        

/*
 
	private void setClusterLabelsMaxIndex() {

        int numNodes = treeModel.getNodeCount();

		int K_int =0;
		for(int i=0; i < numNodes; i++){
			if( (int) indicators.getParameterValue(i) == 1){
				K_int++;
			}
		}
		
        int[] cutNodes = new int[K_int];
 	   int cutNum = 0;
 	   String content = "";
        for(int i=0; i < numNodes; i++){
     	   if( (int) indicators.getParameterValue( i ) ==1 ){
     		   cutNodes[cutNum] = i;
     		   content +=  i +  ",";
     		   cutNum++;
     	   }
     	  
        }
     //   System.out.println(content + " K_int=" + K_int);
        
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
        

        
        
//        int anotherCount = 0;
 //       for(int i=0; i < 20; i++){
  //      	if( ((int) excisionPoints.getParameterValue(i)) ==1)
   //     	anotherCount ++;
  //      }
//        System.out.println();
   //     int maxLabel=0;
    //    for(int i=0; i < numdata; i++){
     //   	if(maxLabel < membership[membershipToClusterLabelIndexes[i]]){
      //  		maxLabel = membership[membershipToClusterLabelIndexes[i]];
       // 	}
       // }
        
       // System.out.println("anotherCount=" + anotherCount + " K_int = " + K_int + " max label=" + maxLabel);
        
        
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


 
 */

/*
  

	private void resetStatusAndBreakpointsGivenCutNodes(int[] testCutNode, int[] onPoints) {
	
		for(int i=0; i < numNodes; i++){
			indicators.setParameterValue(i, 0);
		}
			
	 	int numOn = testCutNode.length;
	 	int countOn=0;
	 	for(int i=0; i < onPoints.length; i++){
			if(countOn < numOn){
				indicators.setParameterValue(onPoints[countOn], 1);
				countOn++;
	 		}
	 		
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
    
 */
        
        
        //
        
        //THE OLD MULTI-STEP
        /*
		else if(operatorSelect == 2){
//			else if(chooseOperator > 0.3){
				

				
				//new move... move the node up and down, with the hope that it promotes mixing.
			 		
		  		int isOn = 0;
		  		int I_selected = -1;
		  		while(isOn == 0){
		  			I_selected = (int) (Math.floor(Math.random()*binSize));
		  			isOn = (int) status.getParameterValue(I_selected);  			
		  		}    	  		
		  		
		 		
		 		
		 		
					//Determining the number of steps from the original site
					int []numStepsFromOrigin = new int[numNodes];
					for(int i=0; i < numNodes; i++){
						numStepsFromOrigin[i] = 100000;
					}
		 		
		  		//System.out.println("Excision point = " + excisionPoints.getParameterValue(I_selected));

		  		
		  			int curElementNumber =(int) breakPoints.getParameterValue(I_selected);
		  			
		  			//curElementNumber = 700;//testing purpose
		  			//maxNodeLevel= 10000000;
		  			
		  			
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
		  		    
		  		    //int xcount=0;
		  		  //System.out.println("root node " + curElement.getNumber());
	 		    while(visitlist.size() > 0){
	 		    	//xcount++;
				
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
			  				if( breakPoints.getParameterValue(i) == site_test){
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
		  					    
	 		    
	 		    //System.out.println("# visit = " + xcount);
	 		    //System.exit(0);
	 		    
	 		   // for(int i=0; i < possibilities.size(); i++){
	 		   // 	System.out.println(possibilities.get(i));
	 		    //}
	 		    //System.out.println(possibilities.size());
	//System.exit(0);
		  			int numPossibleMoves = possibilities.size();  // this number may be different from the max number because some moves may not be legal.
		 			//System.out.print("  num possible moves = " + numPossibleMoves + "\t[");
		  			//for(int i=0; i < numPossibleMoves; i++){
		  			//	System.out.print( possibilities.get(i) + ",");
		  			//}
		  			//System.out.println("]");
		  			//create a list of possible moves
			 		
		  			
		  			int whichMove = (int) (Math.floor(Math.random()*numPossibleMoves));

		  			
		  			//for(int i=0; i < numPossibleMoves; i++){
		  				//System.out.println(movePossibilities[i]);
		  			//}
		  			

		  			//
		  			int site_add = possibilities.get(whichMove);
		  			breakPoints.setParameterValue(I_selected, site_add);
		  			
		  			
		  			howmanyStepsMultiStep = numStepsFromOrigin[site_add];
		  		   
		  			
		  	//    	System.out.println("Chose " + site_add + "  (steps=" + howmanyStepsMultiStep + ")");  		
		//  			System.out.println("propose from " + curElementNumber + " to " + site_add);
		  		
		  		selectedI = I_selected;
		  		
		  		
		  		//calculate the MH requires me to know the (number of) backward moves...
	  			//System.out.println("curElementNumber=" + curElementNumber);
		  		
	  			curElementNumber = site_add;
	  			rootElementNumber = curElementNumber;
	  			//System.out.println("curElementNumber=" + curElementNumber);
	  			curElement = treeModel.getNode(curElementNumber); 
	  			
	  		    visitlist = new LinkedList<NodeRef>();
	  		    fromlist = new LinkedList<NodeRef>();
	  		    nodeLevel = new LinkedList<Integer>();
	  		    
	  		    possibilities = new LinkedList<Integer>();
	  		    
	  		    dummyNode = null;
	  		    visitlist.add(curElement);
	  		    fromlist.add(dummyNode);
	  		    nodeLevel.add(new Integer(0));
	  		    
	  		    
	  		  //System.out.println("root node " + curElement.getNumber());
			    while(visitlist.size() > 0){
			
	  			if(treeModel.getParent(curElement) != null){
	  				//add parent
			  			NodeRef node= treeModel.getParent(curElement);	  		  			
	  				if(fromlist.getFirst() != node){
	  					if( nodeLevel.getFirst() <= maxNodeLevel){
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
	  					if( nodeLevel.getFirst() <= maxNodeLevel){
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
		  				if( breakPoints.getParameterValue(i) == site_test){
		  					hasBeenAdded=1;
		  					break;
		  				}
		  			}
		  			if(hasBeenAdded==0 || curElement.getNumber() == rootElementNumber ){
		  				//System.out.println("to possibilities: add " + site_test);
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

			    
			   // for(int i=0; i < possibilities.size(); i++){
			   // 	System.out.println(possibilities.get(i));
			    //}
			    //System.out.println(possibilities.size());
	//System.exit(0);
		  		
		  		int newStatenumPossibleMoves = possibilities.size();  // this number may be different from the max number because some moves may not be legal.
	  			//System.out.println("num new states moves = " + newStatenumPossibleMoves);
		  		//System.out.println("#  possibilities = " + numPossibleMoves);
		  		//System.out.println("# new possibilities = " + newStatenumPossibleMoves);

		  		
		  		
		  		
	//System.exit(0);
		  		
		  		logHastingRatio = Math.log( (1/ (double)newStatenumPossibleMoves) / (1/ (double)numPossibleMoves)  );
		  		//System.out.println("log hasting ratio = " + logHastingRatio);	
	  			//System.exit(0);
		  		
		  		
		  		
		  		
		  		
		  		
		  		
	//logHastingRatio = 100000;	  		//for testing only - the move should always be accepted
//		  		System.out.println("treeClusterAlg's MH ratio =" + Math.exp(logHastingRatio));
 }
		  		*/


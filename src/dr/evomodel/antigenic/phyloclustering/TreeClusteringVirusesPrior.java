package dr.evomodel.antigenic.phyloclustering;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.Iterator;
import java.util.LinkedList;

import dr.evolution.datatype.Nucleotides;
import dr.evolution.tree.NodeRef;
import dr.evomodel.antigenic.phyloclustering.misc.obsolete.TiterImporter;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.MatrixParameter;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.math.MathUtils;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLObjectParser;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

/**
 * @author Charles Cheung
 * @author Trevor Bedford
 */


public class TreeClusteringVirusesPrior extends AbstractModelLikelihood {
	

    public static final String TREE_CLUSTER_VIRUSES = "treeClusterViruses";
    
    private Parameter muPrecision;
    private Parameter p_on;
    private double initialK;
    
    private int startBase;
    private int endBase;
    private int numSites;
    
	//Variables
    private Parameter indicators;
    private Parameter probSites;
    private MatrixParameter mu; //mu - means
    private Parameter clusterLabels;     //C
    private Parameter clusterLabelsTreeNode;
    private MatrixParameter virusLocations = null;
    private MatrixParameter virusLocationsTreeNode = null;

    private Parameter muMean = null; //obsolete - to be removed
    boolean hasDrift;
    TreeModel treeModel;
    
    int numdata;
    int []membershipToClusterLabelIndexes = null;        
	private double mostRecentTransformedValue = 0;  //keep a copy of the most recent version of transformFactor, to keep track of whether the transformFactor has changed
	private boolean treeChanged = false; //a flag that becomes true when treeModel changes


    private String[] mutationString;
    private LinkedList<Integer>[] mutationList ;
    private LinkedList<Integer>[] causalList ;
    private int[] causalCount;
    private int[] nonCausalCount;

    private Parameter siteIndicators;
    private String gp_prior;
    private double prob00=0.95;
    private double prob11=0.5;


    public LinkedList<Integer>[] getMutationList(){
    	return(mutationList);
    }
    public LinkedList<Integer>[] getCausalList(){
    	return(causalList);
    }
        
	public TreeClusteringVirusesPrior (TreeModel treeModel_in,  
		 				Parameter indicators_in, 
		 				Parameter clusterLabels_in,
		 				Parameter clusterLabelsTreeNode_in,
		 				MatrixParameter mu_in, 
		 				Boolean hasDrift, 
		 				MatrixParameter virusLocations_in,
		 				MatrixParameter virusLocationsTreeNode_in,
		 				Parameter muPrecision, 
		 				Parameter p_onValue	,
		 				double initialK_in,
		 				Parameter muMean,
		 				Parameter probSites_in,
		 				Parameter siteIndicators_in,
		 				int startBase_in,
		 				int endBase_in,
		 				String gp_prior_in,
		 				double prob00_in,
		 				double prob11_in,
		 				double initial_probSiteValue
						){
	 
		super(TREE_CLUSTER_VIRUSES);	
		System.out.println("loading the constructor for TreeClusterViruses");
		
		if(probSites_in == null){
			System.out.println("Antigenic Clustering only");
		}
		else{
			System.out.println("Antigenic Genotype to Phenotype model");
		}

		this.treeModel= treeModel_in;
		
		if(probSites_in != null){
			 this.startBase = startBase_in;
			 this.endBase = endBase_in;
			 this.gp_prior = gp_prior_in;
			 
	 		System.out.println("gp prior = " + gp_prior);
	 		if(gp_prior.compareTo("generic") == 0 ||gp_prior.compareTo("saturated") == 0||
	 				gp_prior.compareTo("shrinkage") == 0||gp_prior.compareTo("correlated") == 0){
	 		}
	 		else{
	 			System.out.println("Prior is incorrectly specified - choose from [generic/saturated/shrinkage/correlated]");
	 			System.exit(0);
	 		}

			treeMutations(); //this routine also calculates the number of sites and checks for error
			this.probSites = probSites_in;
			this.siteIndicators = siteIndicators_in;
			this.prob00 = prob00_in;
			this.prob11 = prob11_in;
	
		}
		//this.K = K_in; //to be deleted
		this.indicators = indicators_in;

		this.clusterLabels = clusterLabels_in;
		this.clusterLabelsTreeNode = clusterLabelsTreeNode_in;
		this.mu = mu_in;
		
		this.hasDrift=hasDrift;
		this.virusLocations = virusLocations_in;
		this.virusLocationsTreeNode = virusLocationsTreeNode_in;
        
		this.muPrecision = muPrecision;
		this.p_on = p_onValue; // this is shared by the clustering-only   and also the saturated prior in genotype to phenotype
		this.initialK = initialK_in;
		this.muMean = muMean;
				
		System.out.println("sigmaSq = " + 1/muPrecision.getParameterValue(0));
		System.out.println("p_on = " + p_onValue.getParameterValue(0));
		
		//numdata = offsets.getSize();
		 int numNodes = treeModel.getNodeCount();
         numdata = virusLocations_in.getColumnDimension();      
        //initialize clusterLabels
        clusterLabels.setDimension(numdata);
         for (int i = 0; i < numdata; i++) {      	
            clusterLabels.setParameterValue(i, 0);
        }
         addVariable(clusterLabels);
         
         //initialize excision points
          indicators.setDimension(treeModel.getNodeCount());
          for(int i=0; i < treeModel.getNodeCount(); i++){
        	  indicators.setParameterValue(i, 0);
          }
          
          

 		 //initialize with the specified number of initial clusters.
 		 indicators.setParameterValue(  treeModel.getRoot().getNumber() , 1 );
 		 for(int k=0; k< (initialK -1); k++){
 				//sample another one.
 				int sampleNew = 1;
 				while(sampleNew ==1){
 					int rSiteIndex = (int) (Math.floor(Math.random()*numNodes));
 					if( (int) indicators.getParameterValue(rSiteIndex) == 0){
 						//clustering only
 						if(probSites_in == null){
 							//success sampling
 							indicators.setParameterValue(rSiteIndex, 1);
 							sampleNew = 0;
 						}
 						//genotype to phenotype
 						else{
 							//added condition to avoid getting 0 probability:
 							if(mutationList[rSiteIndex] != null){				
 								//success sampling
 								indicators.setParameterValue(rSiteIndex, 1);
 								sampleNew = 0;
 							}
 						}//else
 					}
 				}
 				
 		}
          
         addVariable(indicators);
         
         //genotype to phenotype
 		if(probSites_in != null){
	         //int numSites = 330;  //HARDCODED RIGHT NOW
 			if(gp_prior.compareTo("generic") == 0 ){
 				probSites.setDimension(1);
 				probSites.setParameterValue(0, initial_probSiteValue);
 			}
 			else{
 				probSites.setDimension(numSites); //initialize dimension of probSites
 		         //initialize the probability of each site
 		         //double initial_p = 0.05;
 		         for(int k=0; k < numSites; k++){
 		        	 //probSites.setParameterValue(k, initial_p);
 		        	 probSites.setParameterValue(k, initial_probSiteValue );
 		         }  
 			}
	         addVariable(probSites);
	    
	    
	         //MAYBE ONLY INITIALIZE IT IF IT USES THE SHRINKAGE OR correlated PRIOR?
	         siteIndicators.setDimension(numSites);
	         for(int k=0; k < numSites; k++){
	        	 siteIndicators.setParameterValue(k, 1);
	         }
	         addVariable(siteIndicators);
 		}
         
          
          clusterLabelsTreeNode.setDimension(treeModel.getNodeCount());
          addVariable(clusterLabelsTreeNode);
          
         //initialize mu
         mu.setColumnDimension(2);
         mu.setRowDimension( treeModel.getNodeCount()  );  //have a fixed number, although in reality, only K of them are active

         for(int i=0; i < (treeModel.getNodeCount() ); i++){
        	 double zero=0;
        	 mu.getParameter(i).setValue(0, zero);
        	 mu.getParameter(i).setValue(1, zero);
         }
         
         //adding the pre-clustering step.
        // preClustering();

		 addVariable(virusLocations);
		 
		 virusLocationsTreeNode.setColumnDimension(2);  //mds dimension is 2
		 virusLocationsTreeNode.setRowDimension(numNodes);
		 addVariable(virusLocationsTreeNode);
		 
		 addModel(treeModel);
		 addVariable(mu);
		 
		 addVariable(muPrecision);
		 addVariable(p_on);
	
		  numdata = virusLocations.getColumnDimension();
	 
		   //loadInitialMuLocations();
		   //loadIndicators();
	
		   //int []membershipToClusterLabelIndexes = new int[numdata];
		   //setMembershipToClusterLabelIndexes(); //run once in case the tree changes.
		   //setClusterLabelsParameter();
	 
		  //loadIndicators();
		 //System.out.println("Finished loading the constructor for TreeClusteringVirusesPrior");

		 //genotype to phenotype
		if(probSites_in != null){
		 sampleCausativeStates();
		 
		 /*
		 for(int i=0; i< treeModel.getNodeCount(); i++){
			 System.out.print( (int)indicators.getParameterValue(i) + "\t");
			 System.out.print("node i=" + i +":\t");
			 if(mutationList[i] != null){
		    	Iterator itr = mutationList[i].iterator();
		    	Iterator itr2 = causalList[i].iterator();
		    	while(itr.hasNext()){
		    		int curMutation = ((Integer) itr.next()).intValue();
		    		int curCausalState = ((Integer) itr2.next()).intValue();
	 				System.out.print(curMutation + "(" + curCausalState +")" + " ");    		
		    	}
			 }
			 System.out.println("");
			 
		 }
		 */
    	
		}//genotype to phenotype

		 //System.out.println("exit");
		 //System.exit(0);
		 
		 
	}    
 
 

private void treeMutations() {

 	int numNodes = treeModel.getNodeCount();
    // Get sequences
    String[] sequence = new String[numNodes];
    
 // Universal
	String GENETIC_CODE_TABLES ="KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF";
    int numCodons = -1; //needs to set it, when it looks at the sequence

    for(int curIndex = 0; curIndex < numNodes; curIndex ++){
		String ns =  (String) treeModel.getNodeAttribute( treeModel.getNode(curIndex), "states");
		
		if(endBase == -1){
			endBase = ns.length() -1;
		}
		 
		
		 if( (endBase - startBase) % 3 != 0){
			System.out.println("Nucleotide sequence needs to be triplet to convert to codon - check your startbase and endbase");
			System.exit(0);
		 }
		 if(endBase > (ns.length()-1 ) ){
			 System.out.println("the last base cannot be greater than the length of the nucleotide. Exit now.");
			 System.exit(0);
		 }
		 if(startBase > (ns.length()-1 ) ){
			 System.out.println("the start base cannot be greater than the length of the nucleotide. Exit now.");
			 System.exit(0);
		 }
		 if(startBase > endBase){
			 System.out.println("Start base cannot be greater than the end base");
			 System.exit(0);
		 }
		 
		numSites = (endBase - startBase)/3;
	    numCodons = numSites;
   
		//System.out.println("startbase = " + startBase);
		//System.out.println("endbase = " + endBase);
		//System.out.println("numSites = " + numSites);
	   
		ns = ns.substring(startBase, endBase );
		//ns = ns.substring(3+27, ns.length() - 1);
		//System.out.println(ns);
		//System.exit(0);
		
		//numCodons = ns.length()/3;  // or do I care about only 330?

		//System.out.println(numCodons);
		String codonSequence = "";
		for(int codon=0; codon< numCodons; codon++){
			
			int nuc1 =  Nucleotides.NUCLEOTIDE_STATES[ns.charAt(codon*3)];
			int nuc2 =  Nucleotides.NUCLEOTIDE_STATES[ns.charAt(codon*3+1)];
			int nuc3 =  Nucleotides.NUCLEOTIDE_STATES[ns.charAt(codon*3+2)];
			
			int canonicalState = (nuc1 * 16) + (nuc2 * 4) + nuc3;
			
			codonSequence = codonSequence + GENETIC_CODE_TABLES.charAt(canonicalState);
		}
		//System.out.println(codonSequence);
        sequence[curIndex] = codonSequence;
		
    }

    mutationList = new LinkedList[ numNodes];
    mutationString = new String[treeModel.getNodeCount()];

	NodeRef cNode = treeModel.getRoot();
    LinkedList<NodeRef> visitlist = new LinkedList<NodeRef>();
    
    visitlist.add(cNode);
    
    int countProcessed=0;
    while(visitlist.size() > 0){
    	countProcessed++;
    	//assign value to the current node...
    	if(treeModel.getParent(cNode) == null){  //this means it is a root node
    		//visiting the root
    		//System.out.println(cNode.getNumber() + ":\t" + "root");
    	}
    	else{
    		//visiting
    		//System.out.print(cNode.getNumber() + ":\t");

    		//String listMutations = "\"";
    		mutationString[cNode.getNumber()]  = "\"";
    		String nodeState =  sequence[cNode.getNumber()];
    		String parentState =  sequence[treeModel.getParent(cNode).getNumber()];
    		           
    		int count = 0;
    		for(int i=0; i < numCodons; i++){
    			if(nodeState.charAt(i) != parentState.charAt(i)){
    				count++;
    				if(count>1){
    					//System.out.print(",");
    					mutationString[cNode.getNumber()] =  mutationString[cNode.getNumber()] + ",";
    				}
    				//System.out.print(i+1);
    				mutationString[cNode.getNumber()] =  mutationString[cNode.getNumber()] + (i+1);  //i+1 so mutation starts from 1 - 330
    				
    			      // Make sure the list is initialized before adding to it
    			      if (mutationList[cNode.getNumber()] == null) {
    			    	  mutationList[cNode.getNumber()] = new LinkedList<Integer>();
    			      }
    			      mutationList[cNode.getNumber()].add((i+1));
    				
    			}
    			
    			//store in linked list
    		}
    		//System.out.println("");
    		mutationString[cNode.getNumber()]  = mutationString[cNode.getNumber()]  + "\"";
    	}
    	
		//System.out.println(cNode.getNumber() + "\t" +  treeModel.getNodeAttribute(cNode, "states") );

    	
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

	
}



//new version, with probSites determining the indicators
public double getLogLikelihood() {
	
	double N_nodes = (double) treeModel.getNodeCount();
	
	int K_value = 0; //K_int gets updated
	for(int i=0; i < indicators.getDimension(); i++){
		K_value += (int) indicators.getParameterValue(i);
	}
	//System.out.println("K_value" + K_value);

	double logL = 0;
	
	double muVariance = 1/ muPrecision.getParameterValue(0);
	double p_onValue = p_on.getParameterValue(0);
	
	//sync with the current value of p_on
    //for(int k=0; k < treeModel.getNodeCount(); k++){
   	 //probSites.setParameterValue(k, initial_p);
   	// probSites.setParameterValue(k, p_on.getParameterValue(0) );
    //} 
	
	
	double muMeanParameter = muMean.getParameterValue(0);
	//System.out.println("muMeanParameter = " + muMeanParameter);
	
	for(int i=0; i < ( N_nodes ) ; i++){

		double mu_i0 = mu.getParameter(i).getParameterValue(0);
		double mu_i1 = mu.getParameter(i).getParameterValue(1);
		
		//if( (int) indicators.getParameterValue(i) == 1){   //Commented out because I am not using P(mu_i = 0 | I_i = 0) = 1
			logL -=	0.5*(  (mu_i0 - muMeanParameter )*(mu_i0  - muMeanParameter ) + ( mu_i1 )*( mu_i1 )   )/muVariance;
		//}
		//System.out.println(logL);
	}

	
	//clustering
	if(probSites == null){
		logL -= (K_value ) * ( Math.log(2)  + Math.log(Math.PI)+ Math.log(muVariance)  );
	}
	//genotype to phenotype
	else{
		logL -= (N_nodes ) * ( Math.log(2)  + Math.log(Math.PI)+ Math.log(muVariance)  );
		
		//4 priors:
			//generic:  assumes that each $\pi_j = \pi$, where $\pi \sim Beta(\alpha,\beta)$. The posterior $\pi$ is estimated from the MCMC run.
			//Saturated prior: allows each $\pi_j$ to be different. $\pi_j$ is assumed to follow the hierarchical $Beta(\alpha, \beta)$ distribution, where $\alpha$ and $\beta$ are fixed according to plausible prior belief. The $Beta$ distribution is chosen because it is conjugate to the categorical distribution, so Gibbs sampling can be accomplished for $\pi$ (See details in the Implementation). 
			//Shrinkage Prior: models the belief about epitope and non-epitope sites. Here, we define a latent vector of binary indicators for amino acid positions $\delta = (\delta_1, \dots, \delta_L)$, where $\delta_k= 1$ if position $j$ is an epitope site and 0 otherwise. 		
			//correlated prior: models the correlation among amino acid positions being epitopes or not (Figure \ref{correlationPrior}). This prior is motivated by the knowledge that groups of epitope sites tend to occur at adjacent positions.

		double priorContribution = 0;
		if(gp_prior.compareTo("correlated") == 0 ){
			priorContribution = correlatedPriorComputation();
		}
		else if(gp_prior.compareTo("shrinkage") == 0 ){
			priorContribution = shrinkagePriorComputation();
		}
		else if(gp_prior.compareTo("generic") == 0 ){
			priorContribution = genericPriorComputation();
		}
		//saturated
		else if(gp_prior.compareTo("saturated") == 0 ){
			priorContribution = saturatedPriorComputation();
		}
		else{
			System.out.println("Prior unknown. quit now");
			System.exit(0);
		}
	
		logL += priorContribution;
//System.out.println(logL);
// System.exit(0);
	if(logL != Double.NEGATIVE_INFINITY){
		sampleCausativeStates();
	}

	} //end of genotype to phenotype
	return(logL);
}


private double saturatedPriorComputation() {
	double contribution = 0;
	double N_nodes = (double) treeModel.getNodeCount();
	
	double []probMutationSite = new double[numSites]; 
	for(int i=0; i < numSites; i++){
		probMutationSite[i] = probSites.getParameterValue(i);  
 	}
	
	//OMIT THE ROOT NODE BECAUSE THERE IS NO MUTATION ANYWAY
	//(N_nodes - 1)
	for(int i=0; i < (N_nodes-1) ; i ++){
		double prob_Node_i_On = 1;
		double prob_allMutationsOff = 1;
    	if(mutationList[i] != null){
	    	Iterator itr = mutationList[i].iterator();
	    	while(itr.hasNext()){
	    		int curMutation = ((Integer) itr.next()).intValue();
				prob_allMutationsOff = prob_allMutationsOff * (1 - probMutationSite[curMutation -1] );  //offset of 1 
	    	}
    	}  
		
		prob_Node_i_On = 1 - prob_allMutationsOff;
		if( (int) indicators.getParameterValue(i) == 1){
			contribution += Math.log(prob_Node_i_On);
		}
		else{
			contribution += Math.log(1-prob_Node_i_On);
		}
	}
	return(contribution);
}
private double genericPriorComputation() {
	double contribution = 0;
	double N_nodes = (double) treeModel.getNodeCount();
	
	double probMutationSite = probSites.getParameterValue(0); // a single value.. not a vector
	
	//OMIT THE ROOT NODE BECAUSE THERE IS NO MUTATION ANYWAY
	//(N_nodes - 1)
	for(int i=0; i < (N_nodes-1) ; i ++){
		double prob_Node_i_On = 1;
		double prob_allMutationsOff = 1;
    	if(mutationList[i] != null){
	    	Iterator itr = mutationList[i].iterator();
	    	while(itr.hasNext()){
	    		int curMutation = ((Integer) itr.next()).intValue();
				prob_allMutationsOff = prob_allMutationsOff * (1 - probMutationSite );  //offset of 1 
	    	}
    	}  
		
		prob_Node_i_On = 1 - prob_allMutationsOff;
		if( (int) indicators.getParameterValue(i) == 1){
			contribution += Math.log(prob_Node_i_On);
		}
		else{
			contribution += Math.log(1-prob_Node_i_On);
		}
	}
	return(contribution);
}

private double shrinkagePriorComputation() {
	double contribution = 0;
	double N_nodes = (double) treeModel.getNodeCount();
	
	double []probMutationSite = new double[numSites]; 
	for(int i=0; i < numSites; i++){
		probMutationSite[i] = probSites.getParameterValue(i) * siteIndicators.getParameterValue(i) ;  //with the null-ing out the probSites with siteIndicators
 	}
	
	//OMIT THE ROOT NODE BECAUSE THERE IS NO MUTATION ANYWAY
	//(N_nodes - 1)
	for(int i=0; i < (N_nodes-1) ; i ++){
		double prob_Node_i_On = 1;
		double prob_allMutationsOff = 1;
    	if(mutationList[i] != null){
	    	Iterator itr = mutationList[i].iterator();
	    	while(itr.hasNext()){
	    		int curMutation = ((Integer) itr.next()).intValue();
				prob_allMutationsOff = prob_allMutationsOff * (1 - probMutationSite[curMutation -1] );  //offset of 1 
	    	}
    	}  
		
		prob_Node_i_On = 1 - prob_allMutationsOff;
		if( (int) indicators.getParameterValue(i) == 1){
			contribution += Math.log(prob_Node_i_On);
		}
		else{
			contribution += Math.log(1-prob_Node_i_On);
		}
	}
	
	

	int numSignificantSites = 0;
	for(int i=0; i < numSites; i++){
		if( (int) siteIndicators.getParameterValue(i) ==1){
			numSignificantSites++;
		}
	}
	int numNonSignificantSites = numSites - numSignificantSites;
	//contribution +=  (numSignificantSites*Math.log( 0.45   ) + numNonSignificantSites*Math.log( 0.55 )  );
	//contribution +=  (numSignificantSites*Math.log(0.2) + numNonSignificantSites*Math.log(0.8));
	contribution += numSignificantSites*Math.log(p_on.getParameterValue(0)) + numNonSignificantSites*Math.log(1-p_on.getParameterValue(0));
	//System.out.println("#sig= " + numSignificantSites + " #nonsig=" + numNonSignificantSites);
		
	
	return(contribution);
}
private double correlatedPriorComputation() {
	double contribution = 0;
	double N_nodes = (double) treeModel.getNodeCount();
	
	// p^k (1-p)^(numNodes - k)
	//logL += K_value*Math.log( p_onValue ) + (N_nodes - K_value)*Math.log( 1- p_onValue);
	//int numSites = 330; //now have a private variable 
	double []probMutationSite = new double[numSites]; 
	for(int i=0; i < numSites; i++){
		//probMutationSite[i] = probSites.getParameterValue(i);
		probMutationSite[i] = probSites.getParameterValue(i) * siteIndicators.getParameterValue(i) ;  //with the null-ing out the probSites with siteIndicators
 	}
	
	//OMIT THE ROOT NODE BECAUSE THERE IS NO MUTATION ANYWAY
	//(N_nodes - 1)
	for(int i=0; i < (N_nodes-1) ; i ++){
		double prob_Node_i_On = 1;
		double prob_allMutationsOff = 1;
//		for( each mutation in the node i){
//			prob_allMutationsOff = prob_allMutationsOff * (1 - probSites.getParameterValue(curMutation_node_i)); 
//		}
    	if(mutationList[i] != null){
	    	Iterator itr = mutationList[i].iterator();
	    	while(itr.hasNext()){
	    		int curMutation = ((Integer) itr.next()).intValue();
				prob_allMutationsOff = prob_allMutationsOff * (1 - probMutationSite[curMutation -1] );  //offset of 1 
	    	}
    	}

	    
		
		prob_Node_i_On = 1 - prob_allMutationsOff;
		if( (int) indicators.getParameterValue(i) == 1){
			contribution += Math.log(prob_Node_i_On);
		}
		else{
			contribution += Math.log(1-prob_Node_i_On);
		}
		//System.out.println(logL);
		//System.out.println("node=" + i + " prob=" + prob_Node_i_On);
	}
	
	

//	int numSignificantSites = 0;
//	for(int i=0; i < numSites; i++){
//		if( (int) siteIndicators.getParameterValue(i) ==1){
//			numSignificantSites++;
//		}
//	}
//	int numNonSignificantSites = numSites - numSignificantSites;
//	logL +=  (numSignificantSites*Math.log( 0.45   ) 
//					+ numNonSignificantSites*Math.log( 0.55 )  );
	
	//System.out.println("#sig= " + numSignificantSites + " #nonsig=" + numNonSignificantSites);

	//logL +=  (numSignificantSites*Math.log(0.2) + numNonSignificantSites*Math.log(0.8));
	
	

	//transition matrix:
	
	// p(0->0) = 0.9, p(0->1) = 0.2
	// p(1->0) = 0.6  p(1->1) = 0.4

	contribution += Math.log(0.5); //initial
	int numSignificantSites = 0;
	int num_00=0;
	int num_01=0;
	int num_10=0;
	int num_11=0;
	for(int i=1; i < numSites; i++){
		if( (int) siteIndicators.getParameterValue(i-1) ==0){
			if( (int) siteIndicators.getParameterValue(i) ==0){
				num_00++;
			}
			else{
				num_01++;
			}
		}
		else{
			if( (int) siteIndicators.getParameterValue(i) ==0){
				num_10++;
			}
			else{
				num_11++;
			}
		}
	}
	
	double alpha = prob00; //0 stay as 0
	double beta = prob11; // 1 stay as 1
	//		0		1
	//0 	0.9		0.1
	//1		0.5		0.5
	
	contribution +=  ( num_00*Math.log( alpha ) + num_01*Math.log( 1-alpha )  
			  +num_10*Math.log( 1-beta ) + num_11*Math.log(beta ) );
	
	return(contribution);
}
public int[] getCausalCount(){
	return(causalCount);
}

public int[] getNonCausalCount(){
	return(nonCausalCount);
}

public LinkedList<Integer>[] getMutationsPerNode(){
	return(mutationList);
}

public LinkedList<Integer>[] getCausativeStatesPerNode(){
	return(causalList);
}


//the thing is, if the log likelihood is negative inifinity, this shouldn't be called
// because that means a siteInidicator gets flipped off... so then none of the mutation on the node is on..
//this should produce a 0 likelihood.
//so, no point of realizing the causal state if it is going to be rejected

//as a matter of fact, I think causal state should only be updated upon an acceptance move 
public void sampleCausativeStates(){
	    
	
    causalCount = new int[numSites];  //HARD CODED
    nonCausalCount = new int[numSites]; //HARD CODED
    for(int i=0; i < numSites; i++){
    	causalCount[i] = 0;
    	nonCausalCount[i] = 0;
    }
	
	int N_nodes = (int) treeModel.getNodeCount();
	
	//resample the whole set of causal states
    causalList = new LinkedList[ N_nodes];
	
	for(int curNode=0; curNode < (N_nodes-1) ; curNode ++){
		double prob_Node_i_On = 1;
		double prob_allMutationsOff = 1;
//		for( each mutation in the node i){
//			prob_allMutationsOff = prob_allMutationsOff * (1 - probSites.getParameterValue(curMutation_node_i)); 
//		}
		
		if((int) indicators.getParameterValue(curNode) == 0 ){
	    	if(mutationList[curNode] != null){
	    		causalList[curNode] = new LinkedList<Integer>();
		    	Iterator itr = mutationList[curNode].iterator();
		    	while(itr.hasNext()){
		    		int curMutation = ((Integer) itr.next()).intValue();
		    		causalList[curNode].add(new Integer(0));
		    		nonCausalCount[curMutation -1]++;
		    	}
	    	}
		}
		else{   //if indicator is 1... then need to sample to get the causal indicator (the mutation(s) that give(s) causal state =1
	    	if(mutationList[curNode] != null){
	  
	    		//System.out.println("cur node is " + curNode);
	    		
	    		//count the number of mutations that has nonzero probabilities
		    	//Iterator itr_it = mutationList[curNode].iterator();
		    	//int count = 0;
	    	   	//while(itr_it.hasNext()){
		    	//	int curMutation = ((Integer) itr_it.next()).intValue();
		    	//	if(probSites.getParameterValue(curMutation -1) * siteIndicators.getParameterValue(curMutation -1) >0){
		    	//		count++;
		    	//	}
		    	//}
	    	   	//int numMutations = count;  		
	    		int numMutations = mutationList[curNode].size();

	    		causalList[curNode] = new LinkedList<Integer>();
	    		double[] probM = new double[numMutations];
		    	Iterator itr = mutationList[curNode].iterator();
		    	int count = 0;
		    	while(itr.hasNext()){
		    		int curMutation = ((Integer) itr.next()).intValue();
					//prob_allMutationsOff = prob_allMutationsOff * (1 - probSites.getParameterValue(curMutation));
		    		//probM[count]=	probSites.getParameterValue(curMutation -1);
		    		//if(probSites.getParameterValue(curMutation -1) * siteIndicators.getParameterValue(curMutation -1) >0){
		    		  if(gp_prior.compareTo("generic") == 0 ){
		    			   probM[count]=	probSites.getParameterValue(0);
		    		   }
		    		  else if(gp_prior.compareTo("saturated") == 0 ){
		    			   probM[count]=	probSites.getParameterValue(curMutation -1) ;
		    		   }
		    		  else if(gp_prior.compareTo("correlated") == 0 ||gp_prior.compareTo("shrinkage") == 0   ){
		    			   probM[count]=	probSites.getParameterValue(curMutation -1) * siteIndicators.getParameterValue(curMutation -1);
		    		   }
		    			count++;
		    		//}
		    	}
		    	   	
		    	
		    	//System.out.println("numMutations = " + numMutations);
		    	//generate all possibilities - the binary tuples (I think I actually don't have to realize it..)
		    	//int numMutations = 3; //num mutations
		    	double[] probPossibilities = new double[ (int) Math.pow(2,numMutations)];
		        //double[] probM = {0.05, 0.1, 0.2};
		        
		    	for(int kk=0; kk < Math.pow(2, numMutations); kk++){
		    	    int input = kk;
		    		    
		    	    int[] bits = new int[numMutations];
		    	    for (int i = (numMutations-1); i >= 0; i--) {
		    	        bits[i] = ( (input & (1 << i)) != 0  ) ? 1 : 0;;
		    	    }
		        	probPossibilities[kk] = 1;
		    	    for(int curM=0; curM < numMutations; curM++){
		    	    	if(bits[curM] == 1){
		    	    		probPossibilities[kk] = probPossibilities[kk] * probM[curM];
		    	    	}
		    	    	else{
		    	    		probPossibilities[kk] = probPossibilities[kk] * (1- probM[curM]);
		    	    	}
		    	    }
		    	}
		    	probPossibilities[0] = 0; //zero out the 0,0,0
	    	
		    	//sample from the possibilities.. and then reconstruct back the binary tuple
		        int choice = MathUtils.randomChoicePDF(probPossibilities);
		        //System.out.println("choice is " + choice);

   		        int[] bits = new int[numMutations];
		        for (int i = (numMutations-1); i >= 0; i--) {
		            bits[i] = ( (choice & (1 << i)) != 0  ) ? 1 : 0;;
			        causalList[curNode].add(new Integer(bits[i]));
			        
			        if(bits[i] ==1){
			        	causalCount[mutationList[curNode].get(i).intValue() -1 ]++;
			        }
			        else{
			        	nonCausalCount[mutationList[curNode].get(i).intValue()-1]++;
			        }
		        }
		        //System.out.print(choice + " = " + Arrays.toString(bits) + "\t");
		    	    
		     
		      
		     
		    	//System.exit(0);
		    	
		    	
	    	}
		}
	}
	
	/*
	for(int i=0; i < 330; i++){
		System.out.println(i + "\t" + causalCount[i] + " " + nonCausalCount[i]);
	}
	System.out.println("=====================");
	*/
}
    


/*
public double getLogLikelihood() {
	
	double N_nodes = (double) treeModel.getNodeCount();
	
	int K_value = 0; //K_int gets updated
	for(int i=0; i < indicators.getDimension(); i++){
		K_value += (int) indicators.getParameterValue(i);
	}
	//System.out.println("K_value" + K_value);

	double logL = 0;
	
	double muVariance = 1/ muPrecision.getParameterValue(0);
	double p_onValue = p_on.getParameterValue(0);
	double muMeanParameter = muMean.getParameterValue(0);
	//System.out.println("muMeanParameter = " + muMeanParameter);
		//logL -= (K_value ) * ( Math.log(2)  + Math.log(Math.PI)+ Math.log(muVariance)  );
	logL -= (N_nodes ) * ( Math.log(2)  + Math.log(Math.PI)+ Math.log(muVariance)  );
	
	for(int i=0; i < ( N_nodes ) ; i++){

		double mu_i0 = mu.getParameter(i).getParameterValue(0);
		double mu_i1 = mu.getParameter(i).getParameterValue(1);
		
		//if( (int) indicators.getParameterValue(i) == 1){   //Commented out because I am not using P(mu_i = 0 | I_i = 0) = 1
			logL -=	0.5*(  (mu_i0 - muMeanParameter )*(mu_i0  - muMeanParameter ) + ( mu_i1 )*( mu_i1 )   )/muVariance;
		//}
		//System.out.println(logL);
	}

	// p^k (1-p)^(numNodes - k)
	logL += K_value*Math.log( p_onValue ) + (N_nodes - K_value)*Math.log( 1- p_onValue);
	return(logL);
}
*/
   


private void loadIndicators() {

	FileReader fileReader2;
	try {
		//fileReader2 = new FileReader("/Users/charles/Documents/research/antigenic/GenoPheno/Gabriela/results/initialCondition/H3N2.serumLocs.log");
		//fileReader2 = new FileReader("/Users/charles/Documents/researchData/clustering/output/test25/run64/H3N2_mds.breakpoints.log");
	//	fileReader2 = new FileReader("/Users/charles/Documents/researchData/clustering/output/test25/run79/H3N2_mds.indicators.log");
		fileReader2 = new FileReader("/Users/charles/Documents/researchData/clustering/output/test26/run20-test/H3N2_mds.indicatorsStat-120000.log");
		
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
	    	  
	    	  indicators.setParameterValue(i, Double.parseDouble(datavalue[i+1]));
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





private void loadInitialMuLocations() {

	FileReader fileReader2;
	try {
		//fileReader2 = new FileReader("/Users/charles/Documents/research/antigenic/GenoPheno/Gabriela/results/initialCondition/H3N2.serumLocs.log");
		//fileReader2 = new FileReader("/Users/charles/Documents/researchData/clustering/output/test25/run64/H3N2_mds.mu.log");
		fileReader2 = new FileReader("/Users/charles/Documents/researchData/clustering/output/test25/run79/H3N2_mds.mu.log");
		
	      BufferedReader bReader2 = new BufferedReader( fileReader2);

	      String line = null;

	      //skip to the last line
	      String testLine;
	      while ((testLine = bReader2.readLine()) != null){
	    	  line = testLine;
	      }

	     // System.out.println(line);
	      
	      String datavalue[] = line.split("\t");

	      
	       //   System.out.println(serumLocationsParameter.getParameterCount());
	      for (int i = 0; i < mu.getParameterCount(); i++) {
	    	  
	    	  double dim1 = Double.parseDouble(datavalue[i*2+1]);
	    	  double dim2 = Double.parseDouble(datavalue[i*2+2]);
	    	  mu.getParameter(i).setParameterValue(0, dim1);
	    	  mu.getParameter(i).setParameterValue(1, dim2);
	    	 // System.out.println(datavalue[i*2+1]);
//	    	  System.out.println("mu=" + mu.getParameter(i).getParameterValue(0) +"," + mu.getParameter(i).getParameterValue(1));
	   	  
	      }
	      bReader2.close();
	
	} catch (FileNotFoundException e) {
		e.printStackTrace();
	} catch (IOException e) {
		e.printStackTrace();
	}        

	

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



	private void setClusterLabelsParameter() {
		int K_int = 0;
	     for(int i=0; i < treeModel.getNodeCount(); i++){
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
		
//		System.out.println(content);

		
	for(int childNum=0; childNum < treeModel.getChildCount(curElement); childNum++){
		list.addFirst(treeModel.getChild(curElement,childNum));
	}
	}

	return(membership);
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












//=====================================================================================================================
        
        
	public  int factorial(int n) {
	    int fact = 1; // this  will be the result
	    for (int i = 1; i <= n; i++) {
	        fact *= i;
	    }
	    return fact;
	}
	

    public Model getModel() {
                return this;
            }

    public void makeDirty() {
            }

    public void acceptState() {
        // DO NOTHING
    }

    public void restoreState() {
        // DO NOTHING
    }

    public void storeState() {
        // DO NOTHING
    }

    

    protected void handleModelChangedEvent(Model model, Object object, int index) {
    	if(model == treeModel){
    		//System.out.println("==========Tree model changes!!!!!!!!=====");
    		treeChanged = true;
    	}
        else{
        }
    }
    
    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
    }





    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {
    	   	
    	public final static String EXCISIONPOINTS = "excisionPoints";
    	public final static String CLUSTERLABELS = "clusterLabels";
    	public final static String CLUSTERLABELSTREENODE = "clusterLabelsTreeNode";

    	public final static String  MU = "mu";

    	public final static String OFFSETS = "offsets";
    	public final static String VIRUS_LOCATIONS = "virusLocations";
    	public final static String VIRUS_LOCATIONSTREENODE = "virusLocationsTreeNodes";
    	
    	public final static String INDICATORS = "indicators";
    	public final static String PROBSITES = "probSites";
    	public final static String SITEINDICATORS = "siteIndicators";


        boolean integrate = false;
        
        
     //   public final static String MUVARIANCE = "muVariance";
        public final static String MUPRECISION = "muPrecision";
        public final static String PROBACTIVENODE = "probActiveNode";
        
        public final static String INITIALNUMCLUSTERS = "initialK";
        public final static String MUMEAN = "muMean";
        
        public final static String STARTBASE = "startNucleotide";
        public final static String ENDBASE = "endNucleotide";

        public final static String PROB00 = "prob00";
        public final static String PROB11 = "prob11";
        

        
        public final static String INITIAL_PROBSITE_VALUE = "initialProbSite";
        
        public final static String GP_PRIOR_OPTION = "gp_prior";
        
        public String getParserName() {
            return TREE_CLUSTER_VIRUSES;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        	
        	
        		double initialK = 10;
            	if (xo.hasAttribute(INITIALNUMCLUSTERS)) {
            		initialK = xo.getDoubleAttribute(INITIALNUMCLUSTERS);
            	}

        		double prob00 = 0.95;
            	if (xo.hasAttribute(PROB00)) {
            		prob00 = xo.getDoubleAttribute(PROB00);
            	}
        		double prob11 = 0.5;
            	if (xo.hasAttribute(PROB11)) {
            		prob11 = xo.getDoubleAttribute(PROB11);
            	}

            	double initial_probSiteValue = 0.05;
            	if (xo.hasAttribute(INITIAL_PROBSITE_VALUE)) {
            		initial_probSiteValue = xo.getDoubleAttribute(INITIAL_PROBSITE_VALUE);
            	}            	

            	


            	
            	
        		int startBase = 0;
            	if (xo.hasAttribute(STARTBASE)) {
            		startBase = xo.getIntegerAttribute(STARTBASE) - 1; //minus 1 because index begins at 0
            	}
        		
        		int endBase = -1;
            	if (xo.hasAttribute(ENDBASE)) {
            		endBase = xo.getIntegerAttribute(ENDBASE) -1 ; //minus 1 because index begins at 0
            	}
        		
        		String gp_prior = "";
            	if (xo.hasAttribute(GP_PRIOR_OPTION)) {
            		gp_prior = xo.getStringAttribute(GP_PRIOR_OPTION) ; //minus 1 because index begins at 0
            	}
            	
            	
                TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);

                XMLObject cxo = xo.getChild(CLUSTERLABELS);
                Parameter clusterLabels = (Parameter) cxo.getChild(Parameter.class);

                cxo = xo.getChild(CLUSTERLABELSTREENODE);
                Parameter clusterLabelsTreeNode = (Parameter) cxo.getChild(Parameter.class);

                
                cxo = xo.getChild(MU);
                MatrixParameter mu = (MatrixParameter) cxo.getChild(MatrixParameter.class);
                
                cxo=xo.getChild(VIRUS_LOCATIONS);
                MatrixParameter virusLocations =(MatrixParameter) cxo.getChild(MatrixParameter.class);
                
                cxo = xo.getChild(VIRUS_LOCATIONSTREENODE);
                MatrixParameter virusLocationsTreeNode =(MatrixParameter) cxo.getChild(MatrixParameter.class);
                
                cxo = xo.getChild(INDICATORS);
                Parameter indicators = (Parameter) cxo.getChild(Parameter.class);
  
                cxo = xo.getChild(SITEINDICATORS);
                Parameter siteIndicators = null;
                if(cxo != null){
                  siteIndicators = (Parameter) cxo.getChild(Parameter.class);
                }
                
                cxo = xo.getChild(PROBSITES);
                Parameter probSites = null;
                if(cxo != null){
                	probSites = (Parameter) cxo.getChild(Parameter.class);
                }
                
                cxo = xo.getChild(MUPRECISION);
                Parameter muPrecision = (Parameter) cxo.getChild(Parameter.class);
                
                cxo = xo.getChild(PROBACTIVENODE);
                Parameter probActiveNode = (Parameter) cxo.getChild(Parameter.class);
                
		        boolean hasDrift = false;
		
		        cxo = xo.getChild(MUMEAN);
		        Parameter muMean = (Parameter) cxo.getChild(Parameter.class);
		        
		        return new TreeClusteringVirusesPrior(treeModel, indicators, clusterLabels, clusterLabelsTreeNode, mu, hasDrift, virusLocations, virusLocationsTreeNode, muPrecision, probActiveNode, initialK, muMean, probSites, siteIndicators,
		        		startBase, endBase, gp_prior, prob00, prob11, initial_probSiteValue); 
            }

            //************************************************************************
            // AbstractXMLObjectParser implementation
            //************************************************************************

            public String getParserDescription() {
                return "tree clustering viruses";
            }

            public Class getReturnType() {
                return TreeClusteringVirusesPrior.class;
            }

            public XMLSyntaxRule[] getSyntaxRules() {
                return rules;
            }
            
            
            private final XMLSyntaxRule[] rules = {
                    //AttributeRule.newDoubleRule(MUVARIANCE, true, "the variance of mu"),
                    //AttributeRule.newDoubleRule(PROBACTIVENODE, true, "the prior probability of turning on a node"),
            		AttributeRule.newDoubleRule(INITIALNUMCLUSTERS, true, "the initial number of clusters"),
            		AttributeRule.newDoubleRule(STARTBASE, true, "the start base in the sequence to consider in the genotype to phenotype model"),
            		AttributeRule.newIntegerRule(ENDBASE, true, "the end base in the sequence to consider in the genotype to phenotype model"),
            		AttributeRule.newIntegerRule(STARTBASE, true, "the start base in the sequence to consider in the genotype to phenotype model"),
            		AttributeRule.newDoubleRule(PROB00, true, "correlated prior - the probability of staying at state 0 for adjacent siteIndicator"),
            		AttributeRule.newDoubleRule(PROB11, true, "correlated prior - the probability of staying at state 1 for adjacent siteIndicator"),
            		AttributeRule.newStringRule(GP_PRIOR_OPTION, true, "specifying the prior for probSites (and siteIndicators)"),
            		AttributeRule.newDoubleRule(INITIAL_PROBSITE_VALUE, true, "the initial value of the probSite"),
            		
                    new ElementRule(EXCISIONPOINTS, Parameter.class),
                    new ElementRule(CLUSTERLABELS, Parameter.class),
                    new ElementRule(CLUSTERLABELSTREENODE, Parameter.class),
                    new ElementRule(VIRUS_LOCATIONSTREENODE, MatrixParameter.class),
                    new ElementRule(MU, MatrixParameter.class),
                 //   new ElementRule(OFFSETS, Parameter.class),
                    new ElementRule(VIRUS_LOCATIONS, MatrixParameter.class), 
                    new ElementRule(INDICATORS, Parameter.class),
                    //make it so that it isn't required
                    new ElementRule(SITEINDICATORS, Parameter.class, "the indicator of a site having probability greater than 0 of being associated with antigenic transition", true),
                    new ElementRule(PROBSITES, Parameter.class, "the probability that mutation on a site is associated with antigenic transition", true),
                    new ElementRule(TreeModel.class),
                    new ElementRule(MUPRECISION, Parameter.class),
                    new ElementRule(PROBACTIVENODE, Parameter.class),
                    new ElementRule(MUMEAN, Parameter.class)
            };
            
    };

    String Atribute = null;


	public int getNumSites() {
		return numSites;
	}
	public int getStartBase() {
		return startBase;
	}
	public int getEndBase() {
		// TODO Auto-generated method stub
		return endBase;
	}
        
}












//load initial serum location - load the last line
//OBSOLETE WITH THE NEW indicators	
/*
private void loadBreakpoints() {

		FileReader fileReader2;
		try {
			
			//fileReader2 = new FileReader("/Users/charles/Documents/researchData/clustering/output/test23/run4/H3N2_mds.breakPoints.log");
			
			fileReader2 = new FileReader("/Users/charles/Documents/researchData/clustering/output/test25/run5/H3N2_mds.breakpoints.log");
		      BufferedReader bReader2 = new BufferedReader( fileReader2);

		      String line = null;
	
		      //skip to the last line
		      String testLine;
		      while ((testLine = bReader2.readLine()) != null){
		    	  line = testLine;
		      }

		      System.out.println(line);
		      
		      String datavalue[] = line.split("\t");

		      
		       //   System.out.println(serumLocationsParameter.getParameterCount());
		      for (int i = 0; i < treeModel.getNodeCount(); i++) {
		    	  
		    	  breakPoints.setParameterValue(i, Double.parseDouble(datavalue[i+1]));
		    	 // System.out.println(datavalue[i*2+1]);
		    	  System.out.println("indicators=" + breakPoints.getParameterValue(i));
		   	  
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

*/
 
	
//OBSOLETE WITH INDICATORS..
//load initial serum location - load the last line
/*
private void loadStatus() {

		FileReader fileReader2;
		try {
			//fileReader2 = new FileReader("/Users/charles/Documents/research/antigenic/GenoPheno/Gabriela/results/initialCondition/H3N2.serumLocs.log");
			fileReader2 = new FileReader("/Users/charles/Documents/researchData/clustering/output/test23/run4/H3N2_mds.status.log");
			
		      BufferedReader bReader2 = new BufferedReader( fileReader2);

		      String line = null;
	
		      //skip to the last line
		      String testLine;
		      while ((testLine = bReader2.readLine()) != null){
		    	  line = testLine;
		      }

		      System.out.println(line);
		      
		      String datavalue[] = line.split("\t");

		      
		       //   System.out.println(serumLocationsParameter.getParameterCount());
		      for (int i = 0; i < binSize; i++) {
		    	  
		    	  indicators.setParameterValue(i, Double.parseDouble(datavalue[i+1]));
		    	 // System.out.println(datavalue[i*2+1]);
		    	  System.out.println("excisionPoints=" + indicators.getParameterValue(i));
		   	  
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
*/



/*







  
private void preClustering() {
	
	
	int numNodes = treeModel.getNodeCount(); //need to re-evaluate
	
	int []isOccupied = new int[numNodes];
	 //initial bag:
	   //int numNonZeroIndicators = 0;
	   //nonZeroIndexes = new int[numNodes]; //another variable for quick access of the indexes that are turned on - only the first numNonZeroIndicators are meaningful
	   
	   //assumption: numOn << than the numNodes
	   for(int k=0; k< treeModel.getNodeCount(); k++){
	   	//sample another one.
	   	 int sampleNew = 1;
	   	 while(sampleNew ==1){
	   		int rSiteIndex = (int) Math.floor( Math.random()*numNodes );
	   		if(isOccupied[rSiteIndex] == 0){
	   			//success sampling
	   			//System.out.println("rSiteIndex "+  rSiteIndex);
	   			indicators.setParameterValue( rSiteIndex, 1);
	   			isOccupied[rSiteIndex] = 1;
	   			//nonZeroIndexes[numNonZeroIndicators] = rSiteIndex;
	   			//numNonZeroIndicators++; // be filled to be equal to the numOn
	   			sampleNew = 0;
	   		}
	   	 }
	   }
	   

	   int K_int = (int) K.getParameterValue(0);
	   for(int k=0; k< (K_int-1); k++){
		   	//sample another one.
		   	 int sampleNew = 1;
		   	 while(sampleNew ==1){
		   		int rSiteIndex = (int) Math.floor( Math.random()*numNodes );
		   		if(indicators.getParameterValue(rSiteIndex) == 0){
		   			//success sampling
		   			indicators.setParameterValue(rSiteIndex , 1);
		   			//nonZeroIndexes[numNonZeroIndicators] = rSiteIndex;
		   			//numNonZeroIndicators++; // be filled to be equal to the numOn
		   			sampleNew = 0;
		   		}

		   	 }
		 }
		   
	   //for(int i=0; i < binSize; i++){
		//   System.out.println("excision point = " + excisionPoints.getParameterValue(i));
	   //}

	   //for(int i=0; i < numNonZeroIndicators; i++){
	   	//System.out.println(nonZeroIndexes[i]);	
	   //}

	
//       NodeRef node = treeModel.getRoot();
       
      // if(treeModel.isExternal(node)){
    //	   System.out.println("External node");
     //  }
      // else{
    //	   System.out.println("Internal node");
     //  }
	
       
       int[] cutNodes = new int[K_int];
	   int cutNum = 0;
       for(int i=0; i < numNodes; i++){
    	   if( (int) indicators.getParameterValue( i ) ==1 ){
    		   cutNodes[cutNum] =  i;
    		   cutNum++;
    	   }
    	  
       }
         
   //    for(int i=0; i < K_int; i++){
   // 	   System.out.println(cutNodes[i]);
    //   }
       
       int []membership = determine_membership(treeModel, cutNodes, K_int-1);                           
       
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
	   //which is also unnecessary
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
    	   
    	   clusterLabels.setParameterValue( i, membership[membershipToClusterLabelIndexes[i]]);
    //	   System.out.println(vloc.getParameterName() + " i="+ i + " membership=" + (int) clusterLabels.getParameterValue(i));
    	   
    	 //  Parameter v = virusLocations.getParameter(i);
    	  // System.out.println(v.getParameterName());
       }
       
       
   //    System.out.println("Exit now");
   //    System.exit(0);
       

	int numViruses = offsets.getSize();
	System.out.println("# offsets = " + offsets.getSize());
	//for(int i=0; i < offsets.getSize(); i++){
		//System.out.println(offsets.getParameterValue(i));
	//}
	
	
	//Need a routine to convert the membership back to clusterlabels for external nodes..

	
	int maxLabel=0;
	for(int i=0;i< numdata; i++){
		if(maxLabel < (int) clusterLabels.getParameterValue(i)){
			maxLabel = (int) clusterLabels.getParameterValue(i);
		}
	}
	
	
	//now, change the mu..
	for(int i=0; i <= maxLabel; i++){
		//System.out.println(meanYear[i]*beta);
		//mu.getParameter(i).setParameterValue(0, meanYear[i]*beta);//now separate out mu from virusLocation
		mu.getParameter(i).setParameterValue(0, 0);
		mu.getParameter(i).setParameterValue(1, 0);
	}	
	
	
	//System.exit(0);

*/
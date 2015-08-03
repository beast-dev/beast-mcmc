package dr.evomodel.antigenic.phyloClustering;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.ListIterator;
import java.util.Set;
import java.util.logging.Logger;

import dr.evolution.tree.Tree;
import dr.evolution.tree.NodeRef;
import dr.evomodel.antigenic.phyloClustering.misc.obsolete.TiterImporter;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.treelikelihood.TreeTraitParserUtilities;
import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.CompoundParameter;
import dr.inference.model.MatrixParameter;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.inference.operators.MCMCOperator;
import dr.math.GammaFunction;
import dr.math.distributions.MultivariateNormalDistribution;
import dr.math.matrixAlgebra.SymmetricMatrix;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.ElementRule;
import dr.xml.StringAttributeRule;
import dr.xml.XMLObject;
import dr.xml.XMLObjectParser;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

/**
 * @author Charles Cheung
 * @author Trevor Bedford
 */


public class TreeClusteringVirusesPrior extends AbstractModelLikelihood {
	

    public static final String TREE_CLUSTER_VIRUSES = "TreeClusterViruses";
    
    
    private Parameter muPrecision;
    private Parameter p_on;
    private double initialK;
    
	//Variables
    private Parameter indicators;
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
		 				Parameter muMean
						){
	 
		super(TREE_CLUSTER_VIRUSES);	
		System.out.println("loading the constructor for TreeClusterVIruses");

		this.treeModel= treeModel_in;
		//this.K = K_in; //to be deleted
		this.indicators = indicators_in;
		this.clusterLabels = clusterLabels_in;
		this.clusterLabelsTreeNode = clusterLabelsTreeNode_in;
		this.mu = mu_in;
		
		this.hasDrift=hasDrift;
		this.virusLocations = virusLocations_in;
		this.virusLocationsTreeNode = virusLocationsTreeNode_in;
        
		this.muPrecision = muPrecision;
		this.p_on = p_onValue;
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
 						//success sampling
 						indicators.setParameterValue(rSiteIndex, 1);
 						sampleNew = 0;
 					}
 				}
 				
 		}
          
          
          addVariable(indicators);
          
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
		 System.out.println("Finished loading the constructor for ClusterViruses");

	}    
 
 

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

        boolean integrate = false;
        
        
     //   public final static String MUVARIANCE = "muVariance";
        public final static String MUPRECISION = "muPrecision";
        public final static String PROBACTIVENODE = "probActiveNode";
        
        public final static String INITIALNUMCLUSTERS = "initialK";
        public final static String MUMEAN = "muMean";
        
        public String getParserName() {
            return TREE_CLUSTER_VIRUSES;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        	
        	
        		double initialK = 10;
            	if (xo.hasAttribute(INITIALNUMCLUSTERS)) {
            		initialK = xo.getDoubleAttribute(INITIALNUMCLUSTERS);
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
                
                
                cxo = xo.getChild(MUPRECISION);
                Parameter muPrecision = (Parameter) cxo.getChild(Parameter.class);
                
                cxo = xo.getChild(PROBACTIVENODE);
                Parameter probActiveNode = (Parameter) cxo.getChild(Parameter.class);
                
		        boolean hasDrift = false;
		
		        cxo = xo.getChild(MUMEAN);
		        Parameter muMean = (Parameter) cxo.getChild(Parameter.class);
		        
		        return new TreeClusteringVirusesPrior(treeModel, indicators, clusterLabels, clusterLabelsTreeNode, mu, hasDrift, virusLocations, virusLocationsTreeNode, muPrecision, probActiveNode, initialK, muMean); 
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
                    new ElementRule(EXCISIONPOINTS, Parameter.class),
                    new ElementRule(CLUSTERLABELS, Parameter.class),
                    new ElementRule(CLUSTERLABELSTREENODE, Parameter.class),
                    new ElementRule(VIRUS_LOCATIONSTREENODE, MatrixParameter.class),
                    new ElementRule(MU, MatrixParameter.class),
                 //   new ElementRule(OFFSETS, Parameter.class),
                    new ElementRule(VIRUS_LOCATIONS, MatrixParameter.class), 
                    new ElementRule(INDICATORS, Parameter.class),
                    new ElementRule(TreeModel.class),
                    new ElementRule(MUPRECISION, Parameter.class),
                    new ElementRule(PROBACTIVENODE, Parameter.class),
                    new ElementRule(MUMEAN, Parameter.class)
            };
            
    };

    String Atribute = null;
        
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
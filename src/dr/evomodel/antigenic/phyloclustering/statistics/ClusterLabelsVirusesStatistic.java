
package dr.evomodel.antigenic.phyloclustering.statistics;

import java.util.LinkedList;

import dr.evolution.tree.NodeRef;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.*;
import dr.xml.*;

/**
 *  @author Charles Cheung
 * @author Trevor Bedford
 */


public class ClusterLabelsVirusesStatistic extends Statistic.Abstract implements VariableListener {

	private  int max_dim ;
	
	private LinkedList<Double> path  = new LinkedList<Double>();
//	private int[] fromMembership;
	private int[] clusterLabels;
//	private Parameter clusterLabels;
	private Parameter indicators;
	private MatrixParameter virusLocations;
    private TreeModel treeModel;

    public static final String CLUSTERLABELS_STATISTIC = "clusterLabelsVirusesStatistic";
    
    
    int []membershipToClusterLabelIndexes = null;        


    public ClusterLabelsVirusesStatistic(TreeModel tree, Parameter indicators, int max_dim_in, MatrixParameter virusLocations_in) {
    //public PathStatistic(Parameter clusterLabels, TreeModel tree, Parameter indicators, int max_dim_in) {
      //  this.clusterLabels = clusterLabels;
        this.treeModel = tree;
        this.indicators = indicators;
        this.max_dim = max_dim_in;
        this.virusLocations = virusLocations_in;
        
      //  clusterLabels.addParameterListener(this);
        indicators.addParameterListener(this);
        virusLocations.addParameterListener(this);
        
		setMembershipToClusterLabelIndexes(); // if the tree doesn't change, then I really don't have to do it all the time
    }
    


    public int getDimension() {
		int numdata = virusLocations.getColumnDimension();
        return numdata;
    }



    //assume print in order... so before printing the first number, 
    //determine all the nodes that are active.
    public double getStatisticValue(int dim) {

    	if(dim ==0){
    		//Note: if the tree doesn't change, then I don't have to run setMembershipToClusterLabelIndexes every time
    		setClusterLabelsUsingIndicators();
    	}      
       return ((double) clusterLabels[dim]);

    }

    
    
    

	private void setMembershipToClusterLabelIndexes(){

  	   //I suspect this is an expensive operation, so I don't want to do it many times,
  	   //which is also unnecessary  - MAY have to update whenever a different tree is used.
		int numdata = virusLocations.getColumnDimension();
		int numNodes = treeModel.getNodeCount();
         membershipToClusterLabelIndexes = new int[numdata]; 
         clusterLabels = new int[numdata];
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
	
	
	private void setClusterLabelsUsingIndicators(){

        int []membership = determine_membership_v2(treeModel);
		int numdata = virusLocations.getColumnDimension();
        for(int i=0; i < numdata; i++){   
        	clusterLabels[i] = membership[membershipToClusterLabelIndexes[i]] ;
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

	



    
    
    //traverse down the tree, top down, do calculation
    int[] determine_from_membership_v2(TreeModel treeModel){
	    	//note: I set MAX_DIM as the most I would print, but in order to avoid bug, I 
    	//declare the number of nodes as the most active nodes I can have.
    	int[] fromMembership = new int[treeModel.getNodeCount()];
    	for(int i=0; i < treeModel.getNodeCount(); i++){
    		fromMembership[i ] = -1;
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
	    	 if( (int) indicators.getParameterValue(curElement.getNumber() ) == 1) {
	    //		 System.out.print("indicator # " + curElement.getNumber()  + " ");
	    		numClusters++ ;
	    		membership[ curElement.getNumber() ] = numClusters - 1; 
	    		fromMembership[numClusters -1] = membership[ treeModel.getParent(curElement).getNumber()];
	    //		System.out.println("    membership " + (numClusters-1) + " assigned from " + membership[ treeModel.getParent(curElement).getNumber()] );
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
	
	     return(fromMembership);
   }

	
    
    
    
    

	//private LinkedList<Double> setPath() {
	
		    
		 //return(0);
		
	//}



    
    
    
    
    
    
    
    
    
    
    
    public String getDimensionName(int dim) {
    	
		Parameter v = virusLocations.getParameter(dim);
		String curName = v.getParameterName();
    	
        return curName;
    }

    public void variableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        // do nothing
    	//System.out.println("hi got printed");
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

       // public final static String CLUSTERLABELS = "clusterLabels";
        public final static String INDICATORS = "indicators";        
        public final static String MAXDIMSTR = "maxDim";
        public final static String VIRUSLOCATIONS = "virusLocations";        



        public String getParserName() {
            return CLUSTERLABELS_STATISTIC;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            //Parameter clusterLabels = (Parameter) xo.getElementFirstChild(CLUSTERLABELS);
            TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
            Parameter indicators = (Parameter) xo.getElementFirstChild(INDICATORS);
            MatrixParameter virusLocations = (MatrixParameter) xo.getElementFirstChild(VIRUSLOCATIONS);

        	int maxDim = 30;
        	if(xo.hasAttribute(MAXDIMSTR)){
        		maxDim = xo.getIntegerAttribute(MAXDIMSTR);
        	}

           // return new PathStatistic(clusterLabels, treeModel, indicators, maxDim);
        	 return new ClusterLabelsVirusesStatistic( treeModel, indicators, maxDim, virusLocations);

        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element returns a statistic that shifts a matrix of locations by location drift in the first dimension.";
        }

        public Class getReturnType() {
            return ClusterLabelsVirusesStatistic.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            //new ElementRule(CLUSTERLABELS, Parameter.class),
            new ElementRule(TreeModel.class),
            new ElementRule(INDICATORS, Parameter.class),
            new ElementRule(VIRUSLOCATIONS, Parameter.class),
            AttributeRule.newDoubleRule(MAXDIMSTR, true, "the variance of mu"),
        };
    };

    

}


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

public class PathStatistic extends Statistic.Abstract implements VariableListener {

	private  int max_dim ;
	
	private LinkedList<Double> path  = new LinkedList<Double>();
	private int[] fromMembership;
//	private Parameter clusterLabels;
	private Parameter indicators;
    private TreeModel treeModel;

    public static final String PATH_STATISTIC = "pathStatistic";

    public PathStatistic(TreeModel tree, Parameter indicators, int max_dim_in) {
    //public PathStatistic(Parameter clusterLabels, TreeModel tree, Parameter indicators, int max_dim_in) {
      //  this.clusterLabels = clusterLabels;
        this.treeModel = tree;
        this.indicators = indicators;
        this.max_dim = max_dim_in;
        
      //  clusterLabels.addParameterListener(this);
        indicators.addParameterListener(this);
    }
    


    public int getDimension() {
        return max_dim;
    }



    //assume print in order... so before printing the first number, 
    //determine all the nodes that are active.
    public double getStatisticValue(int dim) {

    	if(dim ==0){
    		//fromMembership = determine_from_membership_v2(treeModel);
    		fromMembership = determine_from_membership_v2();
    	}
    	
    	//need to figure out how many K of them.
      
       return ((double) fromMembership[dim]);

    }

    
    
    int[] determine_from_membership_v3(){
    	
    	int[] printFromCluster = new int[treeModel.getNodeCount()];
    	for(int i=0; i < treeModel.getNodeCount(); i++){
    		printFromCluster[i ] = -1;
    	}
    	
    	
 	    int[] fromMembership = new int[treeModel.getNodeCount()];
    	for(int i=0; i < treeModel.getNodeCount(); i++){
    		fromMembership[i ] = -99;
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
 	    fromMembership[0] = -1;
 	          
 	    while(!list.isEmpty()){
 	    	//do things with the current object
 	    	NodeRef curElement = list.pop();
 	    	
 	    	//cluster assignment:
 	    	if(!treeModel.isRoot(curElement)){
 	    		if( (int) indicators.getParameterValue(curElement.getNumber() ) == 1) {
 	    			numClusters++ ;
 	    			membership[ curElement.getNumber() ] = numClusters - 1;
 	    			fromMembership[ curElement.getNumber() ] = membership[treeModel.getParent(curElement).getNumber()];
 	      	 	}
 	    		else{
 	    			//inherit from parent's cluster assignment
 	    			membership[curElement.getNumber()] = membership[treeModel.getParent(curElement).getNumber()]; 
 	    			fromMembership[curElement.getNumber()] = fromMembership[treeModel.getParent(curElement).getNumber()];
 	    		}        	
 	    	}//is not Root
 	    	else{
 	    		fromMembership[ curElement.getNumber()] = -1;
 	    	}
 	
 	    	
 	        for(int childNum=0; childNum < treeModel.getChildCount(curElement); childNum++){
 	        	list.addFirst(treeModel.getChild(curElement,childNum));
 	        }
 	    }
 	
 	    
 	    for(int i=0; i < treeModel.getNodeCount(); i++){
    		if( (int) indicators.getParameterValue(i ) == 1) {
    			printFromCluster[ membership[i] ] = fromMembership[i];
    		}
 	    }
 	    
 	     return(printFromCluster);
	    
    }
    
    
    
    //traverse down the tree, top down, do calculation
    int[] determine_from_membership_v2(){
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
    	String name = "path" + (dim);
        return name;
    }

    public void variableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        // do nothing
    	//System.out.println("hi got printed");
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

       // public final static String CLUSTERLABELS = "clusterLabels";
        public final static String INDICATORS = "indicators";        
        public final static String MAXDIMSTR = "maxDim";


        public String getParserName() {
            return PATH_STATISTIC;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            //Parameter clusterLabels = (Parameter) xo.getElementFirstChild(CLUSTERLABELS);
            TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
            Parameter indicators = (Parameter) xo.getElementFirstChild(INDICATORS);

        	int maxDim = 30;
        	if(xo.hasAttribute(MAXDIMSTR)){
        		maxDim = xo.getIntegerAttribute(MAXDIMSTR);
        	}

           // return new PathStatistic(clusterLabels, treeModel, indicators, maxDim);
        	 return new PathStatistic( treeModel, indicators, maxDim);

        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element returns a statistic that shifts a matrix of locations by location drift in the first dimension.";
        }

        public Class getReturnType() {
            return PathStatistic.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            //new ElementRule(CLUSTERLABELS, Parameter.class),
            new ElementRule(TreeModel.class),
            new ElementRule(INDICATORS, Parameter.class),
            AttributeRule.newDoubleRule(MAXDIMSTR, true, "the variance of mu"),
        };
    };

    

}

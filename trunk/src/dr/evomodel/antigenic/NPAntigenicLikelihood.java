package dr.evomodel.antigenic;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.logging.Logger;

import com.sun.corba.se.impl.orbutil.graph.Node;

import dr.evolution.tree.Tree;
import dr.evolution.tree.NodeRef;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.treelikelihood.TreeTraitParserUtilities;
import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.inference.operators.MCMCOperator;
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

public class NPAntigenicLikelihood extends AbstractModelLikelihood {
	
	  public static final String NP_ANTIGENIC_LIKELIHOOD = "NPAntigenicLikelihood";

	
	
	public NPAntigenicLikelihood (TreeModel treeModel,CompoundParameter traitParameter , Parameter assignments, Parameter clusterVar, Parameter priorMean, Parameter priorVar, double transformFactor ){
		  super(NP_ANTIGENIC_LIKELIHOOD);

		  this.assignments = assignments;
		  this.clusterVar = clusterVar;
		  this.priorVar = priorVar;
		  this.priorMean = priorMean;
	      this.treeModel= treeModel;
		  this.traitParameter= traitParameter;
		  this.transformFactor=transformFactor;
		  
		  this.alpha= 1.0;
		  
	      
	    int dim = traitParameter.getParameter(0).getSize();
	    numdata = traitParameter.getParameterCount();
		   
	     double Data[][] =new double[numdata][dim];
	    for (int i=0; i<numdata; i++){
	    	for (int j=0; j<dim; j++){
	    	    	Data[i][j]= traitParameter.getParameter(i).getParameterValue(j);
	    	}
	    }
	    
	    this.data=Data;
	    
	    depMatrix=new double[numdata][numdata];
	    List<NodeRef> childList = new ArrayList<NodeRef>(); 
	    this.allTips=Tree.Utils.getExternalNodes(treeModel,treeModel.getRoot());
		
	    recursion(treeModel.getRoot(),childList);
	    logCorrectMatrix(transformFactor);
	    printInformtion(depMatrix);
		logDepMatrix =  new double[numdata][numdata];
	    for(int i=0;i<numdata;i++){
	    	for(int j=0;j<i;j++){
	    		logDepMatrix[i][j]=Math.log(depMatrix[i][j]);
	    		logDepMatrix[j][i]=logDepMatrix[j][i];
	    		
	    	}
	    }
	    
	//double[][] depMatrix1 = getMatrixFromTree(transformFactor);
	//printInformtion(depMatrix1);
	 	 

	 for (int i=0; i<numdata; i++){
		 assignments.setParameterValue(i, i);
	 }
	 
	
	       
	       this.logLikelihoodsVector = new double[assignments.getDimension()+1];
	       
	      double[][] var = new double[2][2];
	     var[0][0]= clusterVar.getParameterValue(0)+priorVar.getParameterValue(0);
	       var[1][1]= clusterVar.getParameterValue(0)+priorVar.getParameterValue(0);
	       var[1][0]=0.0;
	       var[0][1] = 0.0;
	      
	      double[][] precision = new SymmetricMatrix(var).inverse().toComponents();
	       
	       double[] mean = new double[2];
	       mean[0]= priorMean.getParameterValue(0);
	       mean[1]= priorMean.getParameterValue(1);
	       
    	   MultivariateNormalDistribution MVN = new MultivariateNormalDistribution(mean, precision);
	       
    	   for(int i=0;i<logLikelihoodsVector.length-1;i++){
	    	   double[] d = new double[2];
	    	   d[0] = data[i][0];
	    	   d[1] = data[i][1];

	    	   logLikelihoodsVector[i]= MVN.logPdf(d);   
	       }
	     
		
	}

	
	
	
	
	
	  public Model getModel() {
	        return this;
	    }

	  
	  
	  public double[] getLogLikelihoodsVector(){
		  return logLikelihoodsVector;
	  }
	  public double[][] getData(){
		  return data;
	  }
	  
	  public double[][] getDepMatrix(){
		  return depMatrix;
	  }
	  

	  public double[][] getLogDepMatrix(){
		  return logDepMatrix;
	  }
	  
	  public Parameter getPriorMean(){
		  return priorMean;
	  }
	  public Parameter getPriorVar(){
		  return priorVar;
	  }
	  public Parameter getClusterVar(){
		  return clusterVar;
	  }
	  
	  public void setLogLikelihoodsVector(int pos, double value){
		  logLikelihoodsVector[pos]=value;
	  }
	  
	 

	  
	  public double getLogLikelihood() {
		  double logL = 0.0;
		  for (int j=0 ; j<assignments.getDimension();j++){
			  if(logLikelihoodsVector[j]!=0){
				  logL +=logLikelihoodsVector[j];
			  }
		  }
		  
		  for (int j=0 ; j<assignments.getDimension();j++){
		if(assignments.getParameterValue(j)==j){
			logL += Math.log(alpha);
		}
		else{logL += Math.log(depMatrix[j][(int) assignments.getParameterValue(j)]);
			
		}
		
		double sumDist=0.0;
		for (int i=0;i<numdata;i++){
			if(i!=j){sumDist += depMatrix[i][j];
			}
			}
		
		  logL-= Math.log(alpha+sumDist);
		  }
		
	
		  
		 
		 return logL;
	  }
	 
	  
	  
	  /* Marc's suggestion on recursion for getting matrix from tree*/
	  
	  void recursion( NodeRef node, List childList){
		 
		  
		 List<NodeRef> leftChildTipList = new ArrayList<NodeRef>(); 
		 List<NodeRef> rightChildTipList = new ArrayList<NodeRef>(); 
		  
		 if(!treeModel.isExternal(node)){
			 recursion(treeModel.getChild(node, 0),leftChildTipList);
			 recursion(treeModel.getChild(node, 1),rightChildTipList);
			 
			
			 double lBranch = treeModel.getBranchLength(treeModel.getChild(node, 0));
			 double rBranch = treeModel.getBranchLength(treeModel.getChild(node, 1));
			 

			 Set<NodeRef> notLeftChildList = new HashSet<NodeRef>();
			 notLeftChildList.addAll(allTips);
			 	for (NodeRef i :leftChildTipList){
			 			notLeftChildList.remove(i);			 
}
			 	Set<NodeRef> notRightChildList = new HashSet<NodeRef>();
				 notRightChildList.addAll(allTips);
			 	for (NodeRef i :rightChildTipList){
			 		notRightChildList.remove(i);			 
}
			 
			 for (NodeRef lChild : leftChildTipList){
				  for (NodeRef Child : notLeftChildList){
				depMatrix[Child.getNumber()][lChild.getNumber()] += lBranch;	 
				depMatrix[lChild.getNumber()][Child.getNumber()] += lBranch;	 
				 }
			 }


			 for (NodeRef rChild : rightChildTipList){
					 for (NodeRef Child : notRightChildList){
					depMatrix[Child.getNumber()][rChild.getNumber()] += rBranch;	 
					depMatrix[rChild.getNumber()][Child.getNumber()] += rBranch;	 
					 }
				 }
			 
			 
			 childList.addAll(leftChildTipList);
			 childList.addAll(rightChildTipList);
		 }
		 else{
			 childList.add(node);
		 }
		 
	  }
	  
	  
	  
	  void logCorrectMatrix(double p){
		  for (int i=0; i<numdata; i++){
		    	for (int j=0; j<i; j++){
		    		depMatrix[i][j]=1/Math.pow(depMatrix[i][j],p);	    
		    		depMatrix[j][i]=depMatrix[i][j];	    
				    }}
		    
	  }
	  
	  
	  
	  
	  // Slow method for computing matrix from tree
	  
	  
	  
	  public double[][] getMatrixFromTree(double p){
		  double[][] Mat = new double[numdata][numdata];
		  
	  for (int i = 0 ; i<numdata; i++){
		  for (int j =0 ; j<i; j++){
			  Mat[i][j] = -p*Math.log(getTreeDist(i,j));			
			  Mat[j][i] = Mat[i][j];
		  }
	  }
	  return Mat;
	  }
	  
	  
	  
	  
	  public double getTreeDist(int i, int j){
		 double dist=0; 
  
		 NodeRef MRCA = findMRCA(i,j);

		 NodeRef Parent = treeModel.getExternalNode(i);
		 while (Parent!=MRCA){
			 dist+=treeModel.getBranchLength(Parent);
			 Parent = treeModel.getParent(Parent);
		 }

		 Parent = treeModel.getExternalNode(j);
		 while (Parent!=MRCA){
			 dist+=treeModel.getBranchLength(Parent);
			 Parent = treeModel.getParent(Parent);
		 }
		  
			 return dist;
	  }
	  
	  
	  private NodeRef findMRCA(int iTip, int jTip) {
	        Set<String> leafNames = new HashSet<String>();
	        leafNames.add(treeModel.getTaxonId(iTip));
	        leafNames.add(treeModel.getTaxonId(jTip));
	        return Tree.Utils.getCommonAncestorNode(treeModel, leafNames);
	    }


	  public void printInformtion(double[][] Mat) {
         StringBuffer sb = new StringBuffer("matrix \n");
         for(int i=0;i <numdata; i++){
        	 sb.append(" \n");
        	 for(int j=0; j<numdata; j++){
        		 sb.append(Mat[i][j]+" \t");
        	 }
         }
         
         Logger.getLogger("dr.evomodel").info(sb.toString()); };
        
        

   	  public void printOrder() {
            StringBuffer sb = new StringBuffer("taxa \n");
            for(int i=0;i <numdata; i++){
           	 sb.append(" \n");
           	 
           		 sb.append(treeModel.getTaxonId(i));
           	             }
            
            Logger.getLogger("dr.evomodel").info(sb.toString()); };
           
         
         
         
         
         public void printInformtion(double x) {
             StringBuffer sb = new StringBuffer("Info \n");
             		 sb.append(x);
             
             Logger.getLogger("dr.evomodel").info(sb.toString()); };
	  
	
	  
	  
	  
	 
	  
	  
	  
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
	        // DO NOTHING
	    }

	    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
	        // DO NOTHING
	    }

	Set<NodeRef> allTips;  
	CompoundParameter traitParameter;  
	double alpha;
	Parameter clusterVar ;
	Parameter priorVar ;
	Parameter priorMean ;
	Parameter assignments;
	TreeModel treeModel;
	String traitName;
	double[][] data;
	double[][] depMatrix;
	double[][] logDepMatrix;
	double[] logLikelihoodsVector;
	int numdata;
	double transformFactor;
	
	
    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {
        public final static String CLUSTER_VAR = "clusterVar";
        public final static String PRIOR_VAR = "priorVar";
        public final static String PRIOR_MEAN = "priorMean";
        public final static String ASSIGNMENTS = "assignments";
        public final static String TRAIT_NAME = "traitName";
        public final static String TRANSFORM_FACTOR = "transformFactor";
        boolean integrate = false;
    	
        
        public String getParserName() {
            return NP_ANTIGENIC_LIKELIHOOD;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
     
        	TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
        	//String traitName = (String) xo.getAttribute(TRAIT_NAME);
        	
	        XMLObject cxo = xo.getChild(CLUSTER_VAR);
	        Parameter clusterVar = (Parameter) cxo.getChild(Parameter.class);

	        cxo = xo.getChild(PRIOR_VAR);
	        Parameter priorVar = (Parameter) cxo.getChild(Parameter.class);

	        cxo = xo.getChild(PRIOR_MEAN);
	        Parameter priorMean = (Parameter) cxo.getChild(Parameter.class);

	        cxo = xo.getChild(ASSIGNMENTS);
	        Parameter assignments = (Parameter) cxo.getChild(Parameter.class);

	       double transformFactor=1.0;
	        if(xo.hasAttribute(TRANSFORM_FACTOR)){
	        	transformFactor = xo.getDoubleAttribute(TRANSFORM_FACTOR);
	        }
	     
	        TreeTraitParserUtilities utilities = new TreeTraitParserUtilities();
            String traitName = TreeTraitParserUtilities.DEFAULT_TRAIT_NAME;
         
            
            TreeTraitParserUtilities.TraitsAndMissingIndices returnValue =
                    utilities.parseTraitsFromTaxonAttributes(xo, traitName, treeModel, integrate);
           // traitName = returnValue.traitName;
            CompoundParameter traitParameter = returnValue.traitParameter;
            
	     
	        
	        
	        
	       
	        
	        return new NPAntigenicLikelihood(treeModel,traitParameter,  assignments, clusterVar, priorMean,priorVar,transformFactor);
	    }

	    //************************************************************************
	    // AbstractXMLObjectParser implementation
	    //************************************************************************

	    public String getParserDescription() {
	        return "conditional likelihood ddCRP";
	    }

	    public Class getReturnType() {
	        return NPAntigenicLikelihood.class;
	    }

	    public XMLSyntaxRule[] getSyntaxRules() {
	        return rules;
	    }

	    private final XMLSyntaxRule[] rules = {
	    		 new StringAttributeRule(TreeTraitParserUtilities.TRAIT_NAME, "The name of the trait for which a likelihood should be calculated"),
	    	
	    		AttributeRule.newDoubleRule(TRANSFORM_FACTOR,true,"p in transformation of distances -p*log(dist)"),
	    		 
	    		 new ElementRule(TreeTraitParserUtilities.TRAIT_PARAMETER, new XMLSyntaxRule[]{
	                        new ElementRule(Parameter.class)
	                }),   
	    		 new ElementRule(PRIOR_VAR,
		    	                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
		           new ElementRule(CLUSTER_VAR,
			                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
   			       new ElementRule(PRIOR_MEAN,
		    	                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
		    	  new ElementRule(ASSIGNMENTS,
				    	                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
				    	                    new ElementRule(TreeModel.class),
				    	                    
	    };
    };


    
  

    String Atribute = null;
	
}

package dr.evomodel.graph.operators;

import java.io.RandomAccessFile;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Random;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.MutableTree.InvalidTreeException;

import dr.evomodel.graph.GraphModel;
import dr.evomodel.graph.Partition;
import dr.evomodel.graph.PartitionModel;
import dr.inference.model.Parameter;
import dr.inference.operators.AbstractCoercableOperator;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.OperatorFailedException;
import dr.math.MathUtils;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

public class GraphModelDimensionOperator extends AbstractCoercableOperator{	

	public static final String GRAPH_MODEL_DIMENSION_OPERATOR = "graphModelDimensionOperator";
	
	private GraphModel graphModel;
	private PartitionModel partitionModel;
	
	private double probBelowRoot = 0.999;    //This gets transformed in the constructor step.
	private double size = 0.0;            //Translates into add probability of 50%
	
	public GraphModelDimensionOperator(CoercionMode mode, GraphModel graphModel,
									   PartitionModel partitionModel) {
		super(mode);
		this.graphModel = graphModel;
		this.partitionModel = partitionModel;
		
		//This is for computational efficiency
        this.probBelowRoot = -Math.log(1 - Math.sqrt(this.probBelowRoot));
		
		super.setWeight(1.0);
	}

	public double doOperation() throws OperatorFailedException {
        double logq = 0;

        try {
        	if (MathUtils.nextDouble() < Double.POSITIVE_INFINITY){//1.0 / (1 + Math.exp(-size))){
            	logq = addOperation() - size;
            }else{
                logq = removeOperation() + size;
            }
        } catch (NoReassortmentEventException nree) {
            throw new OperatorFailedException("");
        }

        assert !Double.isInfinite(logq) && !Double.isNaN(logq);
        
        return logq;
    }

	
	
	private int findPotentialAttachmentPoints(double time, List<NodeRef> list) {
        int count = 0;

        for (int i = 0, n = graphModel.getNodeCount(); i < n; i++) {
            NodeRef nr = graphModel.getNode(i);
            if (!graphModel.isRoot(nr)) {
                if (graphModel.getNodeHeight(nr) < time) {
                    NodeRef left = getParentWrapper(nr, 0);
                    NodeRef right = getParentWrapper(nr, 1);
                    if (graphModel.isBifurcation(nr)) {
                        assert left == right;
                        if (graphModel.getNodeHeight(left) > time) {
                            if (list != null)
                                list.add(nr);
                            count++;
                        }
                    } else {
                        if (graphModel.getNodeHeight(left) > time) {
                            if (list != null)
                                list.add(nr);
                            count++;
                        }
                        if (graphModel.getNodeHeight(right) > time) {
                            if (list != null)
                                list.add(nr);
                            count++;
                        }
                    }
                }
            } else {
                if (graphModel.getNodeHeight(nr) < time) {
                    if (list != null)
                        list.add(nr);
                    count++;
                }
            }
        }
        return count;
    }
	
	private NodeRef getParentWrapper(NodeRef child, int index){
		if(index == 1 && graphModel.getParent(child, index) != null){
				return graphModel.getParent(child, 1); 	
		}
		return graphModel.getParent(child, index);
	}
	
	
	private double addOperation() throws OperatorFailedException{
		System.out.println(graphModel.linkDump());
		
		double logHastings = 0;
        
        //1. Draw some new heights.
        double treeHeight = graphModel.getNodeHeight(graphModel.getRoot());
        double newBifurcationHeight = Double.POSITIVE_INFINITY;
        double newReassortmentHeight = Double.POSITIVE_INFINITY;

        double theta = probBelowRoot / treeHeight;
		
        while (newBifurcationHeight > treeHeight && newReassortmentHeight > treeHeight) {
            newBifurcationHeight = MathUtils.nextExponential(theta);
            newReassortmentHeight = MathUtils.nextExponential(theta);
        }
       
        logHastings += theta * (newBifurcationHeight + newReassortmentHeight) - Math.log(2.0)
                - 2.0 * Math.log(theta) + Math.log(1 - Math.exp(-2.0 * treeHeight * theta));
        
        if (newBifurcationHeight < newReassortmentHeight) {
            double temp = newBifurcationHeight;
            newBifurcationHeight = newReassortmentHeight;
            newReassortmentHeight = temp;
        }
        

        //2. Find the possible re-assortment and bifurcation points.
        ArrayList<NodeRef> potentialBifurcationChildren = new ArrayList<NodeRef>();
        ArrayList<NodeRef> potentialReassortmentChildren = new ArrayList<NodeRef>();

        int totalPotentialReassortmentChildren = findPotentialAttachmentPoints(
                newReassortmentHeight, potentialReassortmentChildren);

        
        assert totalPotentialReassortmentChildren > 0;
        
        //3.  Choose your re-assortment location.
        NodeRef reassortChild = potentialReassortmentChildren.get(MathUtils
                .nextInt(totalPotentialReassortmentChildren));

        NodeRef reassortLeftParent  = getParentWrapper(reassortChild,0);
        NodeRef reassortRightParent = getParentWrapper(reassortChild,1);
        NodeRef reassortSplitParent = reassortLeftParent;
                      
        if(graphModel.isRecombination(reassortChild)){
        	 boolean[] tester = 
        	 {graphModel.getNodeHeight(reassortLeftParent) > newReassortmentHeight,
              graphModel.getNodeHeight(reassortRightParent) > newReassortmentHeight};
        	 
        	 if(tester[0] && tester[1]){
        		 if(MathUtils.nextBoolean()){
        			 reassortSplitParent = reassortRightParent;
        		 }
        	 }else if(tester[1]){
    			 reassortSplitParent = reassortRightParent;
        	 }
        }
        
    	logHastings += Math.log(2.0);      
        
      	graphModel.beginTreeEdit();
        graphModel.removeChild(reassortSplitParent, reassortChild);
		NodeRef newReassortment = graphModel.newNode();
		graphModel.addChild(reassortSplitParent, newReassortment);
		graphModel.addChild(newReassortment, reassortChild);
        
		graphModel.setNodeHeight(newReassortment, newReassortmentHeight);
        
        //4. Choose your bifurcation location.

        int totalPotentialBifurcationChildren = findPotentialAttachmentPoints(
                newBifurcationHeight, potentialBifurcationChildren);

        assert totalPotentialBifurcationChildren > 0;
        
        logHastings += Math.log((double) potentialBifurcationChildren.size() *
                potentialReassortmentChildren.size());

        NodeRef bifurcationChild = potentialBifurcationChildren.get(MathUtils
                .nextInt(potentialBifurcationChildren.size()));


        NodeRef bifurcationLeftParent =  getParentWrapper(bifurcationChild,0);
        NodeRef bifurcationRightParent = getParentWrapper(bifurcationChild,1);
        NodeRef bifurcationSplitParent = bifurcationLeftParent;

        
      if (graphModel.isRecombination(bifurcationChild)) {
    	  boolean[] tester = {graphModel.getNodeHeight(bifurcationLeftParent) > newBifurcationHeight,
    			  graphModel.getNodeHeight(bifurcationRightParent) > newBifurcationHeight};

    	  if (tester[0] && tester[1]) {
    		  if (MathUtils.nextBoolean()) {
     			  bifurcationSplitParent = bifurcationRightParent;
    		  }
    	  } else if (tester[1]){
  			  bifurcationSplitParent = bifurcationRightParent;
    	  }
      }
      
      graphModel.removeChild(bifurcationSplitParent, bifurcationChild);
	  NodeRef newBifurcation = graphModel.newNode();
      graphModel.addChild(bifurcationSplitParent, newBifurcation);
      graphModel.addChild(newBifurcation, newReassortment);
      if(bifurcationChild==newReassortment){
          graphModel.addChild(newBifurcation, newReassortment);
          logHastings -= Math.log(2.0);
      }else{
    	  graphModel.addChild(newBifurcation, bifurcationChild);
      }
      
      graphModel.setNodeHeight(newBifurcation, newBifurcationHeight);
     	
      // pick one of the partitions on newReassortment at random
      // send it up the new recombination edge
      HashSet<Object> parts = ((GraphModel.Node)reassortChild).getObjects(0);
      for(Object o : parts){
          Partition part = (Partition)o;
          boolean b = MathUtils.nextBoolean();
    	  graphModel.addPartition(newReassortment, b ? 0 : 1, part);
		  graphModel.addPartitionUntilCoalescence(newReassortment, part);
      }

      try{
			graphModel.endTreeEdit();
		}catch(InvalidTreeException e){
			
		}

		System.out.println(graphModel.linkDump() + "\n************************");
		
		return 0;
	}
	

	private int findPotentialNodesToRemove(List<NodeRef> list) {
        int count = 0;
        int n = graphModel.getNodeCount();

        for (int i = 0; i < n; i++) {
            NodeRef node = graphModel.getNode(i);
            
            if(!graphModel.isRecombination(node)){
            	if(graphModel.isBifurcation(graphModel.getParent(node,0)) || 
            			graphModel.isBifurcation(graphModel.getParent(node, 1))){
            		if(list != null){
            			list.add(node);
            		}
            		count++;
            	}
            }
        }
        return count;
    }

	
	private double removeOperation() throws OperatorFailedException{
		 double logHastings = 0;
	        
	     // 1. Draw reassortment node uniform randomly

	     ArrayList<NodeRef> potentialNodes = new ArrayList<NodeRef>();
	     int totalPotentials = findPotentialNodesToRemove(potentialNodes);

	     if (totalPotentials == 0){ 
	    	 throw new NoReassortmentEventException();
	     }
	     
	     logHastings += Math.log((double) totalPotentials);

	     NodeRef recNode = potentialNodes.get(MathUtils.nextInt(totalPotentials));
	     return 0;
	}
	     

	private double deleteOperation() throws OperatorFailedException{
		
		NodeRef[] recombinationNodes = graphModel.getNodesByType(GraphModel.NodeType.RECOMBINANT);
		
		NodeRef rNode = recombinationNodes[0];
		NodeRef rNodeParent = graphModel.getParent(rNode);
		NodeRef rNodeGrandParent = graphModel.getParent(rNodeParent);
		NodeRef rChild = graphModel.getChild(rNode, 0);
		
		graphModel.beginTreeEdit();
		
		graphModel.removeChild(rNode, rChild);
		graphModel.removeChild(rNodeGrandParent, rNodeParent);
		graphModel.removeChild(rNodeParent, rNode);
		graphModel.removeChild(rNodeParent, rNode);
		
		
		graphModel.addChild(rNodeGrandParent, rChild);
		
		graphModel.deleteNode(rNode);
		
		graphModel.deleteNode(rNodeParent);
		
		try{
			graphModel.endTreeEdit();
		}catch(InvalidTreeException e){
			
		}
		

		return 0;
	}

	public String getOperatorName() {
		return GRAPH_MODEL_DIMENSION_OPERATOR;
	}

	public double getCoercableParameter() {
		return 0;
	}

	public double getRawParameter() {
		return 0;
	}

	public void setCoercableParameter(double value) {	
		
	}

	public String getPerformanceSuggestion() {
		return "Write better operators";
	}
	
	 public static dr.xml.XMLObjectParser PARSER = new dr.xml.AbstractXMLObjectParser() {

		public String getParserDescription() {
			return null;
		}

		public Class getReturnType() {
			return GraphModelDimensionOperator.class;
		}

		public XMLSyntaxRule[] getSyntaxRules() {
			return null;
		}

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {
			GraphModel graphModel = (GraphModel) xo.getChild(GraphModel.class);
			PartitionModel partitionModel = (PartitionModel) xo.getChild(PartitionModel.class);
			
			final CoercionMode mode = CoercionMode.parseMode(xo);
 			
			return new GraphModelDimensionOperator(mode,graphModel, partitionModel);
		}

		public String getParserName() {
			return GRAPH_MODEL_DIMENSION_OPERATOR;
		}

	 };
	 
	 private class NoReassortmentEventException extends OperatorFailedException {
	        public NoReassortmentEventException(String message) {
	            super(message);
	        }

	        public NoReassortmentEventException() {
	            super("");
	        }

	        private static final long serialVersionUID = 1L;

	  }

}

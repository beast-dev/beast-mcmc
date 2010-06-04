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
	
	private double probBelowRoot = 0.90;    //This gets transformed in the constructor step.
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
        	if (MathUtils.nextDouble() < 1.0 / (1 + Math.exp(-size))){
            	logq = addOperation() - size;
            }else{
                logq = removeOperation() + size;
            }
        } catch (NoReassortmentEventException nree) {
            return Double.NEGATIVE_INFINITY;
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
	    if(delOp==49)
	    	System.out.println("debugme!");
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
        
    	System.out.println(graphModel.linkDump());
    	
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
      System.out.println(graphModel.linkDump());

      
      NodeRef newBifurcation = graphModel.newNode();
      
     
      if(newBifurcationHeight < treeHeight){
    	  graphModel.removeChild(bifurcationSplitParent, bifurcationChild);
    	  graphModel.addChild(bifurcationSplitParent, newBifurcation);
          graphModel.addChild(newBifurcation, newReassortment);
          if(bifurcationChild==newReassortment){
              graphModel.addChild(newBifurcation, newReassortment);
              logHastings -= Math.log(2.0);
          }else{
        	  graphModel.addChild(newBifurcation, bifurcationChild);
          }
          
          graphModel.setNodeHeight(newBifurcation, newBifurcationHeight);
        
      }else{
    	  //add a new root.
    	  
    	  //swap out the root with newBifurcation
    	  NodeRef root = graphModel.getRoot();
    	  NodeRef rootChild1 = graphModel.getChild(root,0);
    	  NodeRef rootChild2 = graphModel.getChild(root,1);
    	  
    	  graphModel.removeChild(root, rootChild1);
    	  graphModel.removeChild(root, rootChild2);
    	  
    	  graphModel.addChild(newBifurcation,rootChild1);
    	  graphModel.addChild(newBifurcation,rootChild2);
    	  
    	  graphModel.setNodeHeight(newBifurcation, graphModel.getNodeHeight(root));
    	  
    	  //link up the root 
    	  graphModel.addChild(root,newBifurcation);    	  
    	  graphModel.addChild(root, newReassortment);    	  
    	  graphModel.setNodeHeight(root, newBifurcationHeight);
      }
      System.out.println(graphModel.linkDump());
      
      	
      // pick one of the partitions on newReassortment at random
      // send it up the new recombination edge
	  int rcp = graphModel.getParent((GraphModel.Node)reassortChild, 0)==newReassortment ? 0 : 1;
      HashSet<Object> parts = ((GraphModel.Node)reassortChild).getObjects(rcp);
      HashSet<Object> insparts = removeExistingPartitions(parts, (GraphModel.Node)newReassortment);
      for(Object o : insparts){
          boolean b = MathUtils.nextBoolean();
    	  graphModel.addPartition(newReassortment, b ? 0 : 1, (Partition)o);
      }
      // propagate partitions up from newBifurcation's children
	  int bcp = graphModel.getParent((GraphModel.Node)bifurcationChild, 0)==newBifurcation ? 0 : 1;
	  if(bifurcationChild==graphModel.getRoot()){bcp=0;}	//hardwire for root node

	  parts = ((GraphModel.Node)bifurcationChild).getObjects(bcp);
      insparts = removeExistingPartitions(parts, (GraphModel.Node)newBifurcation);
      for(Object o : insparts){
    	  graphModel.addPartition(newBifurcation, 0, (Partition)o);
      }
//	  int nrp = graphModel.getParent((GraphModel.Node)newReassortment, 0)==newBifurcation ? 0 : 1;
      parts = ((GraphModel.Node)newReassortment).getObjects(1);
      insparts = removeExistingPartitions(parts, (GraphModel.Node)newBifurcation);
      for(Object o : parts){
    	  graphModel.addPartition(newBifurcation, 0, (Partition)o);
      }
      
      // ensure that no stale partition objects remain 
      clearStalePartitions((GraphModel.Node)reassortSplitParent);
      clearStalePartitions((GraphModel.Node)newReassortment);
      clearStalePartitions((GraphModel.Node)bifurcationSplitParent);

      validateGraph(graphModel);
      try{
			graphModel.endTreeEdit();
		}catch(InvalidTreeException e){
			
		}
	
		return 0.0;
	}

	private void validateGraph(GraphModel graphModel){
		// sanity check that partitions are always associated with a link
	      for(int nI=0; nI < graphModel.getNodeCount(); nI++){
	    	  GraphModel.Node nn = (GraphModel.Node)graphModel.getNode(nI);
	    	  if(graphModel.getParent(nn,0)==null && nn.getObjects(0).size()>0 && graphModel.getRoot()!=nn)
    	    	  System.err.println("Partitions but no parent!!");
	    	  if(graphModel.getParent(nn,1)==null && nn.getObjects(1).size()>0)
    	    	  System.err.println("Partitions but no parent!!");	    		  
	      }		
	      // sanity check that we haven't added too many partitions!
	      for(int nI=0; nI < graphModel.getNodeCount(); nI++){
	    	  GraphModel.Node nn = (GraphModel.Node)graphModel.getNode(nI);
	    	  if(graphModel.getChildCount(nn) > 0 ||
	    			  graphModel.getParent(nn)!=null){
	    	      if(nn.getObjects(0).size() + nn.getObjects(1).size() > 2)
	    	      {
	    	    	  System.err.println("Too many objects at node!!");
	    	      }
	    	  }
	      }
		
	}
	/*
	 * Finds any partitions that already exist in the Node n and removes them from the set of partitions to be added
	 */
	private HashSet<Object> removeExistingPartitions( HashSet<Object> parts, GraphModel.Node n ){
		HashSet<Object> reduced = new HashSet<Object>();
		for(Object o : parts){
			if(!n.getObjects(0).contains(o) && !n.getObjects(1).contains(o))
				reduced.add(o);
		}
		return reduced;
	}
	
	private void clearStalePartitions(GraphModel.Node awalk){
		if(awalk==null)	return;
    	GraphModel.Node c1 = awalk.getChild(0);
    	GraphModel.Node c2 = awalk.getChild(1);
    	HashSet<Object> allobjs = new HashSet<Object>();
    	if(graphModel.getParent(c1, 0)==awalk)	allobjs.addAll(c1.getObjects(0));
    	if(graphModel.getParent(c1, 1)==awalk)	allobjs.addAll(c1.getObjects(1));
    	if(c2!=null){
    		if(graphModel.getParent(c2, 0)==awalk)	allobjs.addAll(c2.getObjects(0));
    		if(graphModel.getParent(c2, 1)==awalk)	allobjs.addAll(c2.getObjects(1));
    	}
    	HashSet<Object> extra1 = new HashSet<Object>();
    	for(Object o : awalk.getObjects(0)){
    		if(!allobjs.contains(o)){
    			extra1.add(o);
    		}
    	}
    	if(extra1.size()>0){
    		for(Object o : extra1)
    			awalk.removeObject(0, o);
    		clearStalePartitions((GraphModel.Node)graphModel.getParent(awalk,0));    		
    	}

    	HashSet<Object> extra2 = new HashSet<Object>();
    	for(Object o : awalk.getObjects(1)){
    		if(!allobjs.contains(o)){
    			extra2.add(o);
    		}
    	}
    	if(extra2.size()>0){
    		for(Object o : extra2)
    			awalk.removeObject(1, o);
    		clearStalePartitions((GraphModel.Node)graphModel.getParent(awalk,1));    		
    	}
	}

	private int findPotentialNodesToRemove(List<NodeRef> list) {
        int count = 0;
        int n = graphModel.getNodeCount();

        for (int i = 0; i < n; i++) {
            NodeRef node = graphModel.getNode(i);
            
            if(graphModel.isRecombination(node)){
            	if(graphModel.isBifurcation(getParentWrapper(node,0)) || 
            			graphModel.isBifurcation(getParentWrapper(node, 1))){
            		if(list != null){
            			list.add(node);
            		}
            		count++;
            	}
            }
        }
        return count;
    }

	static int delOp = 0;
	
	private double removeOperation() throws OperatorFailedException{
		 double logHastings = 0;
	     delOp++;
	     System.out.println("Deletion op: " + delOp);
	     if(delOp==16)
	    	 System.out.println("debugme!");
	     
	     // 1. Draw reassortment node uniform randomly

	     List<NodeRef> potentialNodes = new ArrayList<NodeRef>();
	     int totalPotentials = findPotentialNodesToRemove(potentialNodes);
	     
	     if (totalPotentials == 0){ 
	    	 throw new NoReassortmentEventException();
	     }
	     
	     logHastings += Math.log((double) totalPotentials);

	     NodeRef removeReassortNode = potentialNodes.get(MathUtils.nextInt(totalPotentials));
	     
	     NodeRef removeReassortNode1 = getParentWrapper(removeReassortNode,0);
	     NodeRef removeReassortNode2 = getParentWrapper(removeReassortNode,1);
	     NodeRef removeBifurcationNode = removeReassortNode1;
	     NodeRef keptNode   = removeReassortNode2;
	     NodeRef removeReassortNodeChild = graphModel.getChild(removeReassortNode, 0);
	     
	     
	     
	     if(removeReassortNode1 != removeReassortNode2){
	    	 if(graphModel.isBifurcation(removeReassortNode1) 
	    			 && graphModel.isBifurcation(removeReassortNode2)){
	    		 if(MathUtils.nextBoolean()){
	    			 removeBifurcationNode = removeReassortNode2;
	    			 keptNode = removeReassortNode1;
	    		 }
	    	}else if(!graphModel.isBifurcation(removeReassortNode1)){
	    			 removeBifurcationNode = removeReassortNode2;
	    			 keptNode = removeReassortNode1;
	    	}
	     }
	     
	     
	     NodeRef removeBifurcationNodeParent = getParentWrapper(removeBifurcationNode,0);
	     NodeRef removeBifurcationNodeChild = graphModel.getChild(removeBifurcationNode,0);
	     
	     if(removeBifurcationNodeChild == removeReassortNode){
	    	 removeBifurcationNodeChild = graphModel.getChild(removeBifurcationNode,1);
	     }
	     
	     if(graphModel.isRoot(removeBifurcationNode)){
	    	 return -100;
	     }
	     
	     graphModel.beginTreeEdit();
	     
	     //Unlink model
	     
	     graphModel.removeChild(removeBifurcationNodeParent, removeBifurcationNode);
	     graphModel.removeChild(removeBifurcationNode, removeBifurcationNodeChild);
	     graphModel.removeChild(removeReassortNode,removeReassortNodeChild);
	     graphModel.removeChild(keptNode, removeReassortNode);
	     
	      
	     if(removeReassortNode1 != removeReassortNode2){
	    	 graphModel.addChild(keptNode, removeReassortNodeChild);
	    	 graphModel.addChild(removeBifurcationNodeParent,removeBifurcationNodeChild);
	    	  graphModel.removeChild(removeBifurcationNode,removeReassortNode);
	   	   
	     }else{
	    	 graphModel.addChild(removeBifurcationNodeParent, removeReassortNodeChild);
	     }
	         
	     graphModel.deleteNode(removeReassortNode);
	     graphModel.deleteNode(removeBifurcationNode);
	     
	     clearStalePartitions((GraphModel.Node)removeBifurcationNodeParent);

	     try{
	     graphModel.endTreeEdit();
	     }catch(InvalidTreeException e){
	    	 
	     }
         System.out.println(graphModel.linkDump());
	     
     	validateGraph(graphModel);
	     
	     
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

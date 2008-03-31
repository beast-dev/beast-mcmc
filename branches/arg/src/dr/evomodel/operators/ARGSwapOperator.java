package dr.evomodel.operators;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.logging.Logger;



import dr.evolution.tree.MutableTree;
import dr.evolution.tree.NodeRef;
import dr.evomodel.tree.ARGModel;
import dr.evomodel.tree.ARGModel.Node;
import dr.inference.operators.OperatorFailedException;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.ElementRule;
import dr.xml.StringAttributeRule;
import dr.xml.XMLObject;
import dr.xml.XMLObjectParser;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

/**
 * This method moves the arg model around.  Use of both the 
 * reassortment and bifurcation modes, as well as the event operator,
 * satisfies irreducibility.
 * @author ebloomqu
 *
 */
public class ARGSwapOperator extends SimpleMCMCOperator{
	
	public static final String ARG_SWAP_OPERATOR = "argSwapOperator";
	public static final String SWAP_TYPE = "type";
	public static final String BIFURCATION_SWAP = "bifurcationSwap";
	public static final String REASSORTMENT_SWAP = "reassortmentSwap";
	public static final String DUAL_SWAP = "dualSwap";
	public static final String FULL_SWAP = "fullSwap";
	public static final String NARROW_SWAP = "narrowSwap";
	public static final String DEFAULT = NARROW_SWAP;
	
	private ARGModel arg;
	private String mode;
	
	public ARGSwapOperator(ARGModel arg, String mode, int weight){
		this.arg = arg;
		this.mode = mode;
		
		setWeight(weight);
		
	}
		
	
	public double doOperation() throws OperatorFailedException {
		
		if(mode.equals(NARROW_SWAP)){
			return narrowSwap();
		}
		
		
		if((mode.equals(REASSORTMENT_SWAP) || mode.equals(DUAL_SWAP)) &&
				arg.getReassortmentNodeCount() == 0){
			return 0.0;
		}
				
		ArrayList<NodeRef> bifurcationNodes = new ArrayList<NodeRef>(arg.getNodeCount());
		ArrayList<NodeRef> reassortmentNodes = new ArrayList<NodeRef>(arg.getNodeCount());
					
		setupBifurcationNodes(bifurcationNodes);
		setupReassortmentNodes(reassortmentNodes);
		
		if(mode.equals(BIFURCATION_SWAP)){
			return bifurcationSwap(bifurcationNodes.get(MathUtils.nextInt(bifurcationNodes.size())));
		}else if(mode.equals(REASSORTMENT_SWAP)){
			return reassortmentSwap(reassortmentNodes.get(MathUtils.nextInt(reassortmentNodes.size())));
		}else if(mode.equals(DUAL_SWAP)){
			reassortmentSwap(reassortmentNodes.get(MathUtils.nextInt(reassortmentNodes.size())));
			return bifurcationSwap(bifurcationNodes.get(MathUtils.nextInt(bifurcationNodes.size())));
		}
				
		bifurcationNodes.addAll(reassortmentNodes);
		
		//TODO Change to heapsort method
		Collections.sort(bifurcationNodes,NodeSorter);
		
		for(NodeRef x : bifurcationNodes){
			if(arg.isBifurcation(x)) 
				bifurcationSwap(x);
			else 
				reassortmentSwap(x);
		}
				
		return 0;
	}
	
	private double narrowSwap() throws OperatorFailedException{
		ArrayList<NarrowSwap> possibleSwaps = new ArrayList<NarrowSwap>(arg.getNodeCount());
		findAllNarrowSwaps(possibleSwaps);
		int possibleSwapsBefore = possibleSwaps.size();
		
		assert possibleSwapsBefore > 0;
		
		doNarrowSwap(possibleSwaps.get(MathUtils.nextInt(possibleSwaps.size())));
		
		possibleSwaps.clear();
		findAllNarrowSwaps(possibleSwaps);
						
		return Math.log((double)possibleSwaps.size()/possibleSwapsBefore);
	}
	
	public int findAllNarrowSwaps(ArrayList<NarrowSwap> moves){
    	NodeRef iP = null, j = null, jP = null;
    	ArrayList<NodeRef> nodes = new ArrayList<NodeRef>(arg.getNodeCount());
    	
    	for(int k = 0, n = arg.getNodeCount(); k < n; k++){
    		NodeRef x = arg.getNode(k);
    		if(!arg.isRoot(x) && !arg.isRoot(arg.getParent(x, 0)) 
    				 && !arg.isRoot(arg.getParent(x, 1))){
    			nodes.add(x);
    		}
    	}
    	NarrowSwap a;
    	for(NodeRef i : nodes){
    		for(int k = 0; k < 2; k++){
    			iP = arg.getParent(i,k);
    			for(int m = 0; m < 2; m++){
    				jP = arg.getParent(iP,m);
    				j = arg.getOtherChild(jP, iP);
    				a = new NarrowSwap(i,iP,j,jP);
    				if(validMove(a) && !moves.contains(a)){
    					moves.add(a);
    				}
    			}
    		}
    	}
    	assert moves.size() > 0;
    	return moves.size();
    }
	
	 private boolean validMove(NarrowSwap move){
    	 if (move.j != move.iP && move.i != move.j &&       
                 (arg.getNodeHeight(move.j) < arg.getNodeHeight(move.iP)) && 
                 (arg.getNodeHeight(move.i) < arg.getNodeHeight(move.jP))) {
         	return true;
             
         }
    	 return false;
    }
	
	
	private void doNarrowSwap(NarrowSwap swap) throws OperatorFailedException{
		 arg.beginTreeEdit();

	        boolean iBifurcation = arg.isBifurcation(swap.i);
	        boolean jBifurcation = arg.isBifurcation(swap.j);

	        if (iBifurcation && jBifurcation) {
	            arg.removeChild(swap.iP, swap.i);
	            arg.removeChild(swap.jP, swap.j);
	            arg.addChild(swap.jP, swap.i);
	            arg.addChild(swap.iP, swap.j);
	        } else if (!iBifurcation && !jBifurcation) {
	        	//dont' do anything
	        } else {
	            if (jBifurcation) {
	                NodeRef t = swap.i;
	                NodeRef tP = swap.iP;
	                swap.i = swap.j;
	                swap.iP = swap.jP;
	                swap.j = t;
	                swap.jP = tP;
	            }
	        }

	        try {
	            arg.endTreeEdit();
	        } catch (MutableTree.InvalidTreeException ite) {
	            throw new OperatorFailedException(ite.toString());
	        }
	}
	
	private class NarrowSwap{
    	public NodeRef i;
    	public NodeRef j;
    	public NodeRef iP;
    	public NodeRef jP;
    	
    	public NarrowSwap(NodeRef i, NodeRef iP, NodeRef j, NodeRef jP){
    		this.i = i;
    		this.j = j;
    		this.iP = iP;
    		this.jP = jP;
    		
    	}
    	
    	public boolean equals(Object o){
    		if(!(o instanceof NarrowSwap)) {
				return false;
    		}
    		
    		NarrowSwap move = (NarrowSwap) o;
    		    		
    		if(this.i == move.i && this.j == move.j &&
    				this.iP == move.iP && this.jP == move.jP){
    			return true;
    		}
    		if(this.i == move.j && this.j == move.i &&
    				this.iP == move.jP && this.jP == move.iP){
    			return true;
    		}
    		
    		return false;
    	}
    	
    	public String toString(){
    		return "(" + i.toString() + ", " + iP.toString() + 
    			   ", " + jP.toString() + ", " + j.toString() + ")";
    	}
    	
    }
	
	
	private double bifurcationSwap(NodeRef x){
		Node startNode = (Node) x;
		
		Node keepChild = startNode.leftChild;
		Node moveChild = startNode.rightChild;
		
		if(MathUtils.nextBoolean()){
			keepChild = moveChild;
			moveChild = startNode.leftChild;
		}
				
		ArrayList<NodeRef> possibleNodes = new ArrayList<NodeRef>(arg.getNodeCount());
		
		findNodesAtHeight(possibleNodes,startNode.getHeight());
		
		assert !possibleNodes.contains(startNode);
		assert possibleNodes.size() > 0;
		
		
		Node swapNode = (Node) possibleNodes.get(MathUtils.nextInt(possibleNodes.size()));
		Node swapNodeParent = swapNode.leftParent;
		
		
		arg.beginTreeEdit();
	
		String before = arg.toARGSummary();
		
		
		if(swapNode.bifurcation){
			swapNodeParent = swapNode.leftParent;
			
			arg.singleRemoveChild(startNode, moveChild);
			
			if(swapNodeParent.bifurcation){
				arg.singleRemoveChild(swapNodeParent, swapNode);
				arg.singleAddChild(swapNodeParent, moveChild);
			}else{
				arg.doubleRemoveChild(swapNodeParent, swapNode);
				arg.doubleAddChild(swapNodeParent, moveChild);
			}
			
			arg.singleAddChild(startNode, swapNode);
			
		}else{
			boolean leftSide = true;
			boolean[] sideOk = {swapNode.leftParent.getHeight() > startNode.getHeight(),
								swapNode.rightParent.getHeight() > startNode.getHeight()};
			
			if(sideOk[0] && sideOk[1]){
				if(MathUtils.nextBoolean()){
					swapNodeParent = swapNode.rightParent;
					leftSide = false;
				}
			}else if(sideOk[1]){
				swapNodeParent = swapNode.rightParent;
				leftSide = false;
			}
			
						
			//Double linked parents
			if(swapNode.leftParent == swapNode.rightParent){
				arg.singleRemoveChild(startNode, moveChild);
						
				if(leftSide){
					swapNode.leftParent = null;
					swapNodeParent.leftChild = null;
				}else{
					swapNode.rightParent = null;
					swapNodeParent.rightChild = null;
				}
						
				arg.singleAddChild(startNode, swapNode);
				arg.singleAddChild(swapNodeParent, moveChild);
			}else if(swapNode.leftParent == startNode || swapNode.rightParent == startNode){
				arg.singleRemoveChild(startNode, moveChild);
						
				if(swapNodeParent.bifurcation){
					arg.singleRemoveChild(swapNodeParent,swapNode);
					arg.singleAddChild(swapNodeParent,moveChild);
				}else{
					arg.doubleRemoveChild(swapNodeParent, swapNode);
					arg.doubleAddChild(swapNodeParent, moveChild);
				}
				
				if(startNode.leftChild == null)
					startNode.leftChild = swapNode;
				else
					startNode.rightChild = swapNode;
				
				if(swapNode.leftParent == null)
					swapNode.leftParent = startNode;
				else
					swapNode.rightParent = startNode;
								
			}else{
				arg.singleRemoveChild(startNode, moveChild);
								
				if(swapNodeParent.bifurcation){
					arg.singleRemoveChild(swapNodeParent, swapNode);
					arg.singleAddChild(swapNodeParent, moveChild);
				}else{
					arg.doubleRemoveChild(swapNodeParent, swapNode);
					arg.doubleAddChild(swapNodeParent, moveChild);
				}
				arg.singleAddChild(startNode, swapNode);
			}
		}
		
		arg.pushTreeChangedEvent();
		
		assert nodeCheck();
		
		try{ 
			arg.endTreeEdit(); 
		}catch(MutableTree.InvalidTreeException ite){
			System.out.println(before);
			System.err.println(ite.getMessage());
			System.exit(-1);
		}
		
		return 0;
	}
	
	private double reassortmentSwap(NodeRef x){
		Node startNode = (Node) x;
		Node startChild = startNode.leftChild;
		
		ArrayList<NodeRef> possibleNodes = new ArrayList<NodeRef>(arg.getNodeCount());		
		
		findNodesAtHeight(possibleNodes,startNode.getHeight());
		
		assert !possibleNodes.contains(startNode);
		assert possibleNodes.size() > 0;
		
		Node swapNode = (Node) possibleNodes.get(MathUtils.nextInt(possibleNodes.size()));
		
		Node swapParent = null;
		
		arg.beginTreeEdit();

		if(swapNode.bifurcation){
			swapParent = swapNode.leftParent;
			
			arg.doubleRemoveChild(startNode, startChild);
			
			if(swapParent.bifurcation)
				arg.singleRemoveChild(swapParent, swapNode);
			else
				arg.doubleRemoveChild(swapParent, swapNode);
			
			arg.doubleAddChild(startNode, swapNode);
			
			
			if(startChild.bifurcation){
				startChild.leftParent = swapParent;
				startChild.rightParent = swapParent;
			}else{
				if(startChild.leftParent == null){
					startChild.leftParent = swapParent;
				}else{
					startChild.rightParent = swapParent;
				}
			}
			if(!swapParent.bifurcation){
				swapParent.leftChild = startChild;
				swapParent.rightChild = startChild;
			}else{
				if(swapParent.leftChild == null){
					swapParent.leftChild = startChild;
				}else{
					swapParent.rightChild = startChild;
				}
			}
					
		}else{
			
			boolean leftSide = true;
			boolean[] sideOk = {swapNode.leftParent.getHeight() > startNode.getHeight(),
								swapNode.rightParent.getHeight() > startNode.getHeight()};
			
			swapParent = swapNode.leftParent;
			
			if(sideOk[0] && sideOk[1]){
				if(MathUtils.nextBoolean()){
					leftSide = false;
					swapParent = swapNode.rightParent;
				}
			}else if(sideOk[1]){
				leftSide = false;
				swapParent = swapNode.rightParent;
			}
			
			if(swapNode.leftParent == swapNode.rightParent){
				arg.doubleRemoveChild(startNode, startChild);
				
				if(leftSide){
					swapParent.leftChild = swapNode.leftParent = null;
					swapParent.leftChild = startChild;
					swapNode.leftParent = startNode;
				}else{
					swapParent.rightChild = swapNode.rightParent = null;
					swapParent.rightChild = startChild;
					swapNode.rightParent = startNode;
				}
				
				startNode.leftChild = startNode.rightChild = swapNode;
				
				if(startChild.bifurcation){
					startChild.leftParent = startChild.rightParent = swapParent;
				}else{
					if(startChild.leftParent == null)
						startChild.leftParent = swapParent;
					else
						startChild.rightParent = swapParent;
				}
			}else{
				arg.doubleRemoveChild(startNode, startChild);
				
				if(swapParent.bifurcation)
					arg.singleRemoveChild(swapParent, swapNode);
				else
					arg.doubleRemoveChild(swapParent, swapNode);
				
				startNode.leftChild = startNode.rightChild = swapNode;
				
				if(leftSide)
					swapNode.leftParent = startNode;
				else
					swapNode.rightParent = startNode;
				
				if(swapParent.bifurcation){
					if(swapParent.leftChild == null)
						swapParent.leftChild = startChild;
					else
						swapParent.rightChild = startChild;
				}else{
					swapParent.leftChild = swapParent.rightChild = startChild;
				}
				
				if(startChild.bifurcation){
					startChild.leftParent = startChild.rightParent = swapParent;
				}else{
					if(startChild.leftParent == null)
						startChild.leftParent = swapParent;
					else
						startChild.rightParent = swapParent;
				}
							
			}
			
		}

		
		arg.pushTreeChangedEvent();
		
		try{ 
			arg.endTreeEdit(); 
		}catch(MutableTree.InvalidTreeException ite){
			System.err.println(ite.getMessage());
			System.exit(-1);
		}
		
		
		return 0;
	}
	
	private void setupBifurcationNodes(ArrayList<NodeRef> list){
		for(int i = 0, n = arg.getNodeCount(); i < n; i++){
			NodeRef x = arg.getNode(i);
			if(arg.isInternal(x) && arg.isBifurcation(x) && !arg.isRoot(x)){
				list.add(x);
			}
		}
	}
	
	private void setupReassortmentNodes(ArrayList<NodeRef> list){
		for(int i = 0, n = arg.getNodeCount(); i < n; i++){
			NodeRef x = arg.getNode(i);
			if(arg.isReassortment(x)){
				list.add(x);
			}
		}
	}
	
	
	private void findNodesAtHeight(ArrayList<NodeRef> x, double height){
		for(int i = 0, n = arg.getNodeCount(); i < n; i++){
			Node test = (Node) arg.getNode(i);
			if(test.getHeight() < height){
				if(test.bifurcation){
					if(test.leftParent.getHeight() > height){
						x.add(test);
					}
				}else{
					if(test.leftParent.getHeight() > height){
						x.add(test);
					}
					if(test.rightParent.getHeight() > height){
						x.add(test);
					}
				}
			}
		}
	}
	
	
	public String getOperatorName() {
		return ARG_SWAP_OPERATOR;
	}

	public String getPerformanceSuggestion() {
		return "";
	}
		
	private Comparator<NodeRef> NodeSorter = new Comparator<NodeRef>(){

		public int compare(NodeRef o1, NodeRef o2) {
			double[] heights = {arg.getNodeHeight(o1),arg.getNodeHeight(o2)};
			
			if(heights[0] < heights[1]){
				return -1;
			}else if(heights[0] > heights[1]){
				return 1;
			}
			
			return 0;
		}
	};
	
	public static XMLObjectParser PARSER = new AbstractXMLObjectParser(){

		public String getParserDescription() {
			return "Swaps nodes on a tree";
		}

		public Class getReturnType() {
			return ARGSwapOperator.class;
		}

		private String[] validFormats = {BIFURCATION_SWAP, REASSORTMENT_SWAP,
				DUAL_SWAP, FULL_SWAP, NARROW_SWAP};
		
		private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
			AttributeRule.newIntegerRule(WEIGHT),	
			new StringAttributeRule(SWAP_TYPE,"The mode of the operator",
					validFormats, false),
			new ElementRule(ARGModel.class),
				
		};
		
		public XMLSyntaxRule[] getSyntaxRules() {
			return rules;
		}

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {
			int weight = xo.getIntegerAttribute(WEIGHT);
			
			String mode = DEFAULT;
			if(xo.hasAttribute(SWAP_TYPE)){
				mode = xo.getStringAttribute(SWAP_TYPE);
			}
			
			Logger.getLogger("dr.evomodel").info("Creating ARGSwapOperator: " + mode);
			
			ARGModel arg = (ARGModel) xo.getChild(ARGModel.class);
			return new ARGSwapOperator(arg,mode,weight);
		}

		public String getParserName() {
			return ARG_SWAP_OPERATOR;
		}
		
	};
	
	public boolean nodeCheck() {
		for (int i = 0, n = arg.getNodeCount(); i < n; i++) {
			Node x = (Node) arg.getNode(i);

			if (x.leftParent != x.rightParent &&
					x.leftChild != x.rightChild) {
				return false;
			}
			if (x.leftParent != null) {
				if (x.leftParent.leftChild.getNumber() != i &&
						x.leftParent.rightChild.getNumber() != i)
					return false;
			}
			if (x.rightParent != null) {
				if (x.rightParent.leftChild.getNumber() != i &&
						x.rightParent.rightChild.getNumber() != i)
					return false;
			}
			if (x.leftChild != null) {
				if (x.leftChild.leftParent.getNumber() != i &&
						x.leftChild.rightParent.getNumber() != i)
					return false;
			}
			if (x.rightChild != null) {
				if (x.rightChild.leftParent.getNumber() != i &&
						x.rightChild.rightParent.getNumber() != i)
					return false;
			}
		}

		return true;
	}

}

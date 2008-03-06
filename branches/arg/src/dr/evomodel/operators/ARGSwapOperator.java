package dr.evomodel.operators;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;

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

public class ARGSwapOperator extends SimpleMCMCOperator{
	
	public static final String ARG_SWAP_OPERATOR = "argSwapOperator";
	public static final String SWAP_TYPE = "swapType";
	public static final String BIFURCATION_SWAP = "bifurcationSwap";
	public static final String REASSORTMENT_SWAP = "reassortmentSwap";
	public static final String DUAL_SWAP = "dualSwap";
	public static final String FULL_SWAP = "fullSwap";
	public static final String DEFAULT_MODE = BIFURCATION_SWAP;
	
	private ARGModel arg;
	private String mode;
	
	public ARGSwapOperator(ARGModel arg, String mode, int weight){
		this.arg = arg;
		this.mode = mode;
		
		setWeight(weight);
	}
	
	
	
	public double doOperation() throws OperatorFailedException {
		
		ArrayList<NodeRef> bifurcationNodes = new ArrayList<NodeRef>(arg.getNodeCount());
		ArrayList<NodeRef> reassortmentNodes = new ArrayList<NodeRef>(arg.getNodeCount());
		
		if(mode.equals(BIFURCATION_SWAP)){
			setupBifurcationNodes(bifurcationNodes);
			
			return bifurcationSwap(bifurcationNodes.get(MathUtils.nextInt(bifurcationNodes.size())));
		}else if(mode.equals(REASSORTMENT_SWAP)){
			if(arg.getReassortmentNodeCount() == 0)
				return 0.0;
			setupReassortmentNodes(reassortmentNodes);
			
			return reassortmentSwap(reassortmentNodes.get(MathUtils.nextInt(reassortmentNodes.size())));
		}else if(mode.equals(DUAL_SWAP)){
			if(arg.getReassortmentNodeCount() == 0)
				return 0.0;
			setupBifurcationNodes(bifurcationNodes);
			setupReassortmentNodes(reassortmentNodes);
			
			bifurcationSwap(bifurcationNodes.get(MathUtils.nextInt(bifurcationNodes.size())));
			return reassortmentSwap(reassortmentNodes.get(MathUtils.nextInt(reassortmentNodes.size())));
		}
		
		setupBifurcationNodes(bifurcationNodes);
		setupReassortmentNodes(reassortmentNodes);
		
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
	
	private double bifurcationSwap(NodeRef x){
		Node node = (Node) x;
		
		Node keepChild = node.leftChild;
		Node moveChild = node.rightChild;
		
		if(MathUtils.nextBoolean()){
			keepChild = moveChild;
			moveChild = node.leftChild;
		}
		
		return 0;
	}
	
	private double reassortmentSwap(NodeRef x){
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
	
	public XMLObjectParser PARSER = new AbstractXMLObjectParser(){

		public String getParserDescription() {
			return "Swaps nodes on a tree";
		}

		public Class getReturnType() {
			return ARGSwapOperator.class;
		}

		private String[] validFormats = {BIFURCATION_SWAP, REASSORTMENT_SWAP,
				DUAL_SWAP, FULL_SWAP};
		
		private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
			AttributeRule.newIntegerRule(WEIGHT),	
			new StringAttributeRule(SWAP_TYPE,"The mode of the operator",
					validFormats, true),
			new ElementRule(ARGModel.class),
				
		};
		
		public XMLSyntaxRule[] getSyntaxRules() {
			return rules;
		}

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {
			int weight = xo.getIntegerAttribute(WEIGHT);
			
			String mode = DEFAULT_MODE;
			if(xo.hasAttribute(SWAP_TYPE)){
				mode = xo.getStringAttribute(SWAP_TYPE);
			}
			
			ARGModel arg = (ARGModel) xo.getChild(ARGModel.class);
			return new ARGSwapOperator(arg,mode,weight);
		}

		public String getParserName() {
			return ARG_SWAP_OPERATOR;
		}
		
	};

}

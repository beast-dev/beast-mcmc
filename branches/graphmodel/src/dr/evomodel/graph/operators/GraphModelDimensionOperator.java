package dr.evomodel.graph.operators;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.MutableTree.InvalidTreeException;
import dr.evomodel.graph.GraphModel;
import dr.evomodel.graph.PartitionModel;
import dr.inference.operators.AbstractCoercableOperator;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.OperatorFailedException;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

public class GraphModelDimensionOperator extends AbstractCoercableOperator{	
	private GraphModel graphModel;
	private PartitionModel partitionModel;
	
	public GraphModelDimensionOperator(CoercionMode mode, GraphModel graphModel,
									   PartitionModel partitionModel) {
		super(mode);
		this.graphModel = graphModel;
		this.partitionModel = partitionModel;
		
		super.setWeight(1.0);
	}

	public static final String GRAPH_MODEL_DIMENSION_OPERATOR = "graphModelDimensionOperator";
	
	public double doOperation() throws OperatorFailedException {
		if(graphModel.getNodeCount() > 7){
			return deleteOperation();
		}
		
		return addOperation();
		
		
	}
	
	private double addOperation() throws OperatorFailedException{
		
		
		
		graphModel.beginTreeEdit();
		
		NodeRef leaf = graphModel.getExternalNode(0);
		
		NodeRef leafParent = graphModel.getParent(leaf);
		
		NodeRef newNode1 = graphModel.newNode();
		NodeRef newNode2 = graphModel.newNode();
		
		graphModel.removeChild(leafParent, leaf);
		
		graphModel.addChild(newNode1, leaf);
		graphModel.addChild(newNode2,newNode1);
		graphModel.addChild(newNode2,newNode1);
		graphModel.addChild(leafParent, newNode2);
		
		graphModel.addPartition(newNode1, partitionModel.getSiteRange(0));
		graphModel.addPartition(newNode2, partitionModel.getSiteRange(1));
		
		graphModel.setNodeHeight(newNode1, graphModel.getNodeHeight(leafParent)/3.0);
		graphModel.setNodeHeight(newNode2, graphModel.getNodeHeight(leafParent)/3.0*2.0);
		
		try{
			graphModel.endTreeEdit();
		}catch(InvalidTreeException e){
			
		}
		
		return -10000;
	}
	
	private double deleteOperation() throws OperatorFailedException{
		
		NodeRef[] recombinationNodes = graphModel.getNodesByType(2);
		
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

}

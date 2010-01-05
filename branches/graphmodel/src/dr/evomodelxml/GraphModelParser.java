package dr.evomodelxml;

import java.util.logging.Logger;

import dr.evolution.tree.Tree;
import dr.evomodel.graph.GraphModel;
import dr.evomodel.graph.PartitionModel;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

public class GraphModelParser extends TreeModelParser{

	public GraphModelParser(){
		super();
		
		rules = null; //Write some rules
	}
	
	public String getParserDescription() {
		return "This element represents a model of a graph.";
	}

	public Class getReturnType() {
		return GraphModel.class;
	}

	public XMLSyntaxRule[] getSyntaxRules() {
		return rules;
	}

	public Object parseXMLObject(XMLObject xo) throws XMLParseException {
		Tree tree = (Tree) xo.getChild(Tree.class);
		PartitionModel partitionModel = (PartitionModel) xo.getChild(PartitionModel.class);
		
		GraphModel graphModel = new GraphModel(xo.getId(), tree, partitionModel);
		
		Logger.getLogger("dr.evomodel").info("Creating the graph model, '" + xo.getId() + "'");

		super.replaceParameters(xo, graphModel);
		
		graphModel.setupHeightBounds();
		
		Logger.getLogger("dr.evomodel").info("  initial tree topology = " + Tree.Utils.uniqueNewick(graphModel, graphModel.getRoot()));
        Logger.getLogger("dr.evomodel").info("  tree height = " + graphModel.getNodeHeight(graphModel.getRoot()));
		
		return graphModel;
	}

	public String getParserName() {
		return GraphModel.GRAPH_MODEL;
	}
	

}

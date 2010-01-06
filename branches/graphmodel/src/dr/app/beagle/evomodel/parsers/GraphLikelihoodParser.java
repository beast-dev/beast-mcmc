package dr.app.beagle.evomodel.parsers;

import dr.app.beagle.evomodel.treelikelihood.PartialsRescalingScheme;
import dr.evomodel.graph.GraphModel;
import dr.evomodel.graph.PartitionModel;
import dr.evomodel.graphlikelihood.GraphLikelihood;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Likelihood;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

public class GraphLikelihoodParser extends TreeLikelihoodParser{

	public String getParserName(){
		return GraphLikelihood.GRAPH_LIKELIHOOD;
	}
	
	public Class getReturnType() {
        return GraphLikelihood.class;
    }
	
	public XMLSyntaxRule[] getSyntaxRules(){
		return new XMLSyntaxRule[]{
				
		};
	}
	
	public Likelihood setupLikelihood(boolean useAmbiguities, 
    		int instanceCount, TreeModel treeModel, PartialsRescalingScheme scalingScheme,
    		XMLObject xo) throws XMLParseException{
			
		PartitionModel pm = (PartitionModel) xo.getChild(PartitionModel.class);
		
		return new GraphLikelihood((GraphModel)treeModel, pm, null, useAmbiguities, 
				true, true, false, false);
	}
}

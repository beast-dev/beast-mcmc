package dr.evomodel.graph;

import dr.inference.model.Statistic;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.XMLObject;
import dr.xml.XMLObjectParser;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

public class RecombinationNodeStatistic extends Statistic.Abstract {

	public static final String RECOMBINATION_NODE_STATISTIC = "recombinationNodeCount";
	public GraphModel graphModel;
	
	public RecombinationNodeStatistic(GraphModel graphModel){
		this.graphModel = graphModel;
	}
	
	public int getDimension() {
		return 1;
	}

	public double getStatisticValue(int dim) {
		return (graphModel.getNodesByType(GraphModel.NodeType.RECOMBINANT)).length;
	}
	
	public static XMLObjectParser PARSER = new AbstractXMLObjectParser(){

		public String getParserDescription() {
			return "";
		}

		@SuppressWarnings("unchecked")
		public Class getReturnType() {
			return RecombinationNodeStatistic.class;
		}

		public XMLSyntaxRule[] getSyntaxRules() {
			return null;
		}

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {
			GraphModel gm = (GraphModel) xo.getChild(GraphModel.class);
			return new RecombinationNodeStatistic(gm);
		}

		public String getParserName() {
			return RecombinationNodeStatistic.RECOMBINATION_NODE_STATISTIC;
		}
		
	};
	
	

}

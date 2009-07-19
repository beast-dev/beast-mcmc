package dr.evomodel.arg;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.arg.ARGModel.Node;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Parameter.ChangeType;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.XMLObject;
import dr.xml.XMLObjectParser;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

public class ARGRelaxedClock extends AbstractModel implements BranchRateModel{
	
	public static final String ARG_LOCAL_CLOCK = "argLocalClock";
	public static final String PARTITION = "partition";
	
	private Parameter globalRateParameter;
	
	private ARGModel arg;
	private int partition;
		
	
	public ARGRelaxedClock(String name) {
		super(name);
	}
	
	public ARGRelaxedClock(String name, ARGModel arg, int partition,Parameter rate){
		super(name);
		
		this.arg=arg;
		this.partition = partition;
		
		globalRateParameter = rate;
		
		addModel(arg);
		addParameter(rate);
	}

	protected void acceptState() {
			
	}

	protected void handleModelChangedEvent(Model model, Object object, int index) {
		//do nothing
	}

	
	protected void handleParameterChangedEvent(Parameter parameter, int index,ChangeType type) {
		//do nothing
	}

	
	protected void restoreState() {
		
	}

	
	protected void storeState() {
		
	}

	public double getBranchRate(Tree tree, NodeRef nodeRef) {
	
		Node treeNode = (Node)nodeRef;
		Node argNode = (Node)treeNode.mirrorNode;
		
		
		return globalRateParameter.getParameterValue(0)*argNode.getRate(partition);
	}

	public String getAttributeForBranch(Tree tree, NodeRef node) {
		return Double.toString(getBranchRate(tree, node));
	}

	public String getBranchAttributeLabel() {
		 return "rate";
	}
	
	public static XMLObjectParser PARSER = new AbstractXMLObjectParser(){

		public String getParserDescription() {
			return null;
		}

		public Class getReturnType() {
			
			return ARGRelaxedClock.class;
		}

		public XMLSyntaxRule[] getSyntaxRules() {
			
			return null;
		}

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {

			ARGModel arg = (ARGModel)xo.getChild(ARGModel.class);
			
			int partition = xo.getAttribute(PARTITION, 0);
			
			Parameter rate = (Parameter)xo.getChild(Parameter.class);
			
			return new ARGRelaxedClock("",arg,partition,rate);
		}

		public String getParserName() {
			return ARG_LOCAL_CLOCK;
		}
		
	};

}

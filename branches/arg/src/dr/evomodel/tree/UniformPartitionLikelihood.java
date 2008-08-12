package dr.evomodel.tree;

import java.util.logging.Logger;

import dr.evomodel.operators.NewerARGEventOperator;
import dr.inference.model.Likelihood;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLObjectParser;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

public class UniformPartitionLikelihood extends Likelihood.Abstract{

	public static final String UNIFORM_PARTITION_LIKELIHOOD = "uniformPartitionLikelihood";
	private double logPartitionNumber; //Transformed initially for computational reasons
	private ARGModel arg;
	
	public UniformPartitionLikelihood(ARGModel arg) {
		super(arg);
		this.arg = arg;
		
		if(arg.isRecombinationPartitionType()){
			logPartitionNumber = -Math.log(arg.getNumberOfPartitions() - 1);
		}else{
			logPartitionNumber = -(arg.getNumberOfPartitions() - 1)*NewerARGEventOperator.LOG_TWO;
		}
	}
	
	public double calculateLogLikelihood() {
		return logPartitionNumber*arg.getReassortmentNodeCount();
	}
	
	public double getLogPartitionNumber(){
		return logPartitionNumber;
	}
	
	public static XMLObjectParser PARSER = new AbstractXMLObjectParser(){

		public String getParserDescription() {
			return "Provides a uniform prior for partitions";
		}

		public Class getReturnType() {
			return UniformPartitionLikelihood.class;
		}

		public XMLSyntaxRule[] getSyntaxRules() {
			return new XMLSyntaxRule[]{
				new ElementRule(ARGModel.class),
			};
		}

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {
			Logger.getLogger("dr.evomodel").info("Creating " + UNIFORM_PARTITION_LIKELIHOOD);
			
			return new UniformPartitionLikelihood((ARGModel)xo.getChild(ARGModel.class));
		}

		public String getParserName() {
			return UNIFORM_PARTITION_LIKELIHOOD;
		}
		
	};

}

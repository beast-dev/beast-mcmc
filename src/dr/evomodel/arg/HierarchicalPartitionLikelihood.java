package dr.evomodel.arg;

import dr.inference.distribution.DistributionLikelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Parameter.ChangeType;
import dr.math.MathUtils;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLObjectParser;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

public class HierarchicalPartitionLikelihood extends ARGPartitionLikelihood{

	public static final String HIERARCHICAL_PARTITION_LIKELIHOOD = "hierarchicalPartitionLikelihood";
	
	private Parameter probabilities;
	
	public HierarchicalPartitionLikelihood(String id, ARGModel arg, Parameter probs) {
		super(id, arg);
		
		this.probabilities = probs;
		
		addParameter(probs);
	}

	public double[] generatePartition() {
		double[] partition = new double[getNumberOfPartitionsMinusOne() + 1];
		
		partition[0] = 0.0;
		
		for(int i = 0; i < partition.length; i++)
			partition[i] = 0.0;
		
		while(UniformPartitionLikelihood.arraySum(partition) == 0.0){
			for(int i = 1; i < partition.length; i++){
				if(MathUtils.nextDouble() < probabilities.getParameterValue(i-1)){
					partition[i] = 1.0;
				}else{
					partition[i] = 0.0;
				}
			}
		}
		
		
		return partition;
	}

	public double getLogLikelihood(double[] partition) {
		double logLike = 0;
		
		for(int i = 1; i < partition.length; i++){
			if(partition[i] == 1.0){
				logLike += Math.log(probabilities.getParameterValue(i-1));
			}else{
				logLike += Math.log(1 - probabilities.getParameterValue(i-1));
			}
		}
		
		
		return logLike;
	}

	protected void acceptState() {
		// nothing to do!
	}

	protected void handleModelChangedEvent(Model model, Object object, int index) {
		// has no submodels
	}

	protected void handleParameterChangedEvent(Parameter parameter, int index,
			ChangeType type) {
		// I'm lazy, so I compute after each step :)
	}

	protected void restoreState() {
		//nothing to restore!
	}

	@Override
	protected void storeState() {
		//nothing to store
		
	}
	
	public static XMLObjectParser PARSER = new AbstractXMLObjectParser(){

		public String getParserDescription() {
			return null;
		}

		public Class getReturnType() {
			return PoissonPartitionLikelihood.class;
		}

		public XMLSyntaxRule[] getSyntaxRules() {
			return new XMLSyntaxRule[]{
				AttributeRule.newDoubleRule(DistributionLikelihood.MEAN, false),
				new ElementRule(ARGModel.class,false),
			};
		}

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {
			String id = "";
			if(xo.hasId())
				id = xo.getId();
			
			ARGModel arg = (ARGModel)xo.getChild(ARGModel.class);
		
			Parameter values = (Parameter)xo.getChild(Parameter.class);
			
			if(values.getDimension() != arg.getNumberOfPartitions() - 1){
				throw new XMLParseException("The dimension of the parameter must equal the number of partitions minus 1 ");
			}
			
			if(arg.isRecombinationPartitionType()){
				throw new XMLParseException(ARGModel.TREE_MODEL + " must be of type " + ARGModel.REASSORTMENT_PARTITION);
			}
			
			return new HierarchicalPartitionLikelihood(id,arg,values);
		}

		public String getParserName() {
			return HIERARCHICAL_PARTITION_LIKELIHOOD;
		}
		
	};

}

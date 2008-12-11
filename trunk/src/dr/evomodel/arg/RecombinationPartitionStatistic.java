package dr.evomodel.arg;

import dr.evomodel.arg.operators.ARGAddRemoveEventOperator;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Statistic;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.XMLObject;
import dr.xml.XMLObjectParser;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;


public class RecombinationPartitionStatistic extends Statistic.Abstract{

	public final static String RECOMBINATION_PARTITION_STATISTIC = "partitionStatistic";
	private int dimension;
	private ARGModel arg;
	
	public RecombinationPartitionStatistic(String id, ARGModel arg){
		
		setId(id);
		
		this.arg = arg;
		
		this.dimension = arg.getNumberOfPartitions() - 1;
	}
	
	
	public int getDimension() {
		return dimension;
	}

	public double getStatisticValue(int dim) {
		
		int number = 0;
		
		CompoundParameter x = arg.getPartitioningParameters();
		
		for(int i = 0, n= arg.getReassortmentNodeCount() ; i < n; i++){
			double[] z = x.getParameter(i).getParameterValues();
			
			if(ARGAddRemoveEventOperator.arraySum(z) == (double)(dimension - dim)){
				number++;
			}
		}
		
		return (double)number;
	}
	
	public static XMLObjectParser PARSER = new AbstractXMLObjectParser(){

		public String getParserDescription() {
			return null;
		}

		public Class getReturnType() {
			return RecombinationPartitionStatistic.class;
		}

		public XMLSyntaxRule[] getSyntaxRules() {
			return null;
		}

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {
			String id = xo.getId();
			
			ARGModel arg = (ARGModel)xo.getChild(ARGModel.class);
			
			return new RecombinationPartitionStatistic(id,arg);
		}

		public String getParserName() {
			return RECOMBINATION_PARTITION_STATISTIC;
		}
		
	};

}

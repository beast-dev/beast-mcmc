package dr.evomodel.arg;

import java.util.ArrayList;

import dr.inference.model.Statistic;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLObjectParser;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;


public class ARGDistinctTreeCountStatistic extends Statistic.Abstract{

	public static final String ARG_TREE_COUNT = "argTreeCount";
	public static final String PRINT_TREES = "printTrees";
	
	private ARGModel arg;
	private int numberOfPartitions;
	
	public ARGDistinctTreeCountStatistic(ARGModel arg){
		this.arg = arg;
		numberOfPartitions = arg.getNumberOfPartitions();
	}
	
	public int getDimension() {
		
		return 1;
	}

	
	public double getStatisticValue(int dim){
		ArrayList<String> listOfTrees = new ArrayList<String>(numberOfPartitions/10);
		
		for(int i = 0; i < numberOfPartitions; i++){
			ARGTree tree = new ARGTree(arg,i);
			
			String newick = tree.getNewickNoBranches();
			if(!listOfTrees.contains(newick)){
				listOfTrees.add(newick);
			}
		}
		
		
		return listOfTrees.size();
	}
	
	
	public String getTrees(){
		ArrayList<String> listOfTrees = new ArrayList<String>(numberOfPartitions/10);
		
		for(int i = 0; i < numberOfPartitions; i++){
			ARGTree tree = new ARGTree(arg,i);
			
			String newick = tree.getNewickNoBranches();
			if(!listOfTrees.contains(newick)){
				listOfTrees.add(newick);
			}
		}
		
		
		return listOfTrees.toString();
	}
	
	public static XMLObjectParser PARSER = new AbstractXMLObjectParser(){

		public String getParserDescription() {
		
			return "Provides number of distinct loci trees";
		}

		public Class getReturnType() {
			return ARGDistinctTreeCountStatistic.class;
		}

		public XMLSyntaxRule[] getSyntaxRules() {
			return new XMLSyntaxRule[]{
					new ElementRule(ARGModel.class),
			};
		}
		
		public Object parseXMLObject(XMLObject xo) throws XMLParseException {
			ARGModel arg = (ARGModel) xo.getChild(ARGModel.class);
			
			return new ARGDistinctTreeCountStatistic(arg);
		}

		public String getParserName() {
			return ARG_TREE_COUNT;
		}
		
	};

	
	
}

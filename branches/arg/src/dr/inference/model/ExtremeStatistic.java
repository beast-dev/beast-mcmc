package dr.inference.model;

import dr.xml.AbstractXMLObjectParser;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLObjectParser;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

public class ExtremeStatistic extends Statistic.Abstract{
	public static final String EXTREME_STATISTIC = "extremeStatistic";
	private Statistic a;
	
	public int getDimension() {
		return 2;
	}

	private double[] getExtreme(){
		double[] extreme = {a.getStatisticValue(0),a.getStatisticValue(0)};
		
		for(int i = 1; i < a.getDimension(); i++){
			double c = a.getStatisticValue(i);
			if(extreme[0] < c){
				extreme[0] = c;
			}
			if(extreme[1] > c){
				extreme[1] = c;
			}
			
		}
		return extreme;
	}
	
	public double getStatisticValue(int dim) {
		double[] e = getExtreme();
		if(dim == 0){
			return e[0];
		}
		return e[1];
	}
	
	public ExtremeStatistic(Statistic a){
		this.a = a;
	}
	
	public static XMLObjectParser PARSER = new AbstractXMLObjectParser(){

		public String getParserDescription() {
			return "This statistic finds the maximum and minimum of values" +
					" in a particular statistic.";
		}
		
		public Class getReturnType() {
			return ExtremeStatistic.class;
		}

		public XMLSyntaxRule[] getSyntaxRules() {
			return rules;
		}
		
		private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
			new ElementRule(Statistic.class)
		};

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {
			Statistic a = (Statistic)xo.getChild(Statistic.class);
			return new ExtremeStatistic(a);
		}

		public String getParserName() {
			return EXTREME_STATISTIC;
		}
		
		
	};
	
	

}

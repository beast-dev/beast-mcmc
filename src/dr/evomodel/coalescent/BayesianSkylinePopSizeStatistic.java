package dr.evomodel.coalescent;

import dr.inference.model.Statistic;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLObjectParser;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

public class BayesianSkylinePopSizeStatistic extends Statistic.Abstract{

	public double time;
	public BayesianSkylineLikelihood bsl;
	public static final String TIME ="time";
	public static final String BAYESIAN_SKYLINE_POP_SIZE_STATISTIC = 
		"generalizedSkylinePopSizeStatistic";
	
	public BayesianSkylinePopSizeStatistic(double time, 
			BayesianSkylineLikelihood bsl){
		this.time = time;
		this.bsl = bsl;
	}
	
	public int getDimension() {
		return 1;
	}

	public double getStatisticValue(int dim) {
		double[] heights = bsl.getGroupHeights();
		double[] sizes = bsl.getPopSizeParameter().getParameterValues();
		
		
		for(int i = 0; i < heights.length; i++){
			if(this.time < heights[i]){
				return sizes[i];
			}
		}
		
		return Double.NaN;
	}
	
	public static XMLObjectParser PARSER = new AbstractXMLObjectParser(){

		public String getParserDescription() {
			return "The pop sizes at the given times";
		}

		public Class getReturnType() {
			return GMRFPopSizeStatistic.class;
		}

		public XMLSyntaxRule[] getSyntaxRules() {
			return null;
		}

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {
		
			double time = xo.getDoubleAttribute(TIME);
			
			
			BayesianSkylineLikelihood bsl = 
				(BayesianSkylineLikelihood)xo.getChild(BayesianSkylineLikelihood.class);
						
			return new BayesianSkylinePopSizeStatistic(time, bsl);
		}

		public String getParserName() {
			return BAYESIAN_SKYLINE_POP_SIZE_STATISTIC;
		}
		
	};

}

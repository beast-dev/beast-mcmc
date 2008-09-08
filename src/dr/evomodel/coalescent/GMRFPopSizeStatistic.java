package dr.evomodel.coalescent;

import dr.inference.model.Statistic;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLObjectParser;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

public class GMRFPopSizeStatistic extends Statistic.Abstract{
	
	private GMRFSkyrideLikelihood gsl;
	private double[] time;
	
	public final static String TIMES = "time";
	public final static String GMRF_POP_SIZE_STATISTIC = "gmrfPopSizeStatistic";
	
	public GMRFPopSizeStatistic(double[] time, GMRFSkyrideLikelihood gsl){
		this.gsl = gsl;
		this.time = time;
	}
	
	public int getDimension() {
		return time.length;
	}

	public double getStatisticValue(int dim) {
		double[] coalescentHeights = gsl.getCoalescentIntervalHeights();
		double[] popSizes = gsl.getPopSizeParameter().getParameterValues();
		
		assert popSizes.length == coalescentHeights.length;
		
		for(int i = 0; i < coalescentHeights.length; i++){
			if(time[dim] < coalescentHeights[i]){
				return popSizes[i];
			}
		}
				
		return -99;
	}
	
	public static XMLObjectParser PARSER = new AbstractXMLObjectParser(){

		public String getParserDescription() {
			return "The pop sizes at the given times";
		}

		public Class getReturnType() {
			return GMRFPopSizeStatistic.class;
		}

		public XMLSyntaxRule[] getSyntaxRules() {
			return new XMLSyntaxRule[]{
				AttributeRule.newDoubleArrayRule(TIMES, false),
				new ElementRule(GMRFSkyrideLikelihood.class),
			};
		}

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {
		
			double[] times = xo.getDoubleArrayAttribute(TIMES);
			
			if(times.length == 0){
				throw new XMLParseException(TIMES + " must contain at least one time");
			}
			for(double a : times){
				if(a < 0){
					throw new XMLParseException(a + " is invalid time. It must be greater than 0");
				}
			}
			
			GMRFSkyrideLikelihood gsl = (GMRFSkyrideLikelihood)xo.getChild(GMRFSkyrideLikelihood.class);
						
			return new GMRFPopSizeStatistic(times, gsl);
		}

		public String getParserName() {
			return GMRF_POP_SIZE_STATISTIC;
		}
		
	};

}

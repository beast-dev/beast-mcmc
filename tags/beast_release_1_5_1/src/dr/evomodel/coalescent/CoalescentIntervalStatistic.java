package dr.evomodel.coalescent;

import dr.inference.model.Statistic;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLObjectParser;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;


public class CoalescentIntervalStatistic extends Statistic.Abstract{

	public static final String COALESCENT_INTERVAL_STATISTIC = "coalescentIntervalStatistic";
	
	private GMRFSkyrideLikelihood acl;
	private int dimension;
	
	public CoalescentIntervalStatistic(GMRFSkyrideLikelihood acl){
		this.acl = acl;
		
		dimension = acl.getCoalescentIntervalHeights().length;
	}
	
	public int getDimension() {
		return dimension;
	}

	public double getStatisticValue(int dim) {
			return acl.getCoalescentIntervalHeights()[dim];
	}

	public static XMLObjectParser PARSER = new AbstractXMLObjectParser(){

		public String getParserDescription() {
			return "";
		}

		public Class getReturnType() {
			return CoalescentIntervalStatistic.class;
		}

		public XMLSyntaxRule[] getSyntaxRules() {
			return new XMLSyntaxRule[]{
					new ElementRule(GMRFSkyrideLikelihood.class),
			};
		}

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {
			return new CoalescentIntervalStatistic(
					(GMRFSkyrideLikelihood)xo.getChild(GMRFSkyrideLikelihood.class));
		}

		public String getParserName() {
			return COALESCENT_INTERVAL_STATISTIC;
		}
		
	};
	

}

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
	
	private OldAbstractCoalescentLikelihood acl;
	
	public CoalescentIntervalStatistic(OldAbstractCoalescentLikelihood acl){
		this.acl = acl;
	}
	
	public int getDimension() {
		return acl.getIntervalCount();
	}

	public double getStatisticValue(int dim) {
			return acl.getInterval(dim);
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
					new ElementRule(OldAbstractCoalescentLikelihood.class),
			};
		}

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {
			return new CoalescentIntervalStatistic(
					(OldAbstractCoalescentLikelihood)xo.getChild(OldAbstractCoalescentLikelihood.class));
		}

		public String getParserName() {
			return COALESCENT_INTERVAL_STATISTIC;
		}
		
	};
	

}

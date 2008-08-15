package dr.evomodel.coalescent;

import dr.inference.model.Statistic;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLObjectParser;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

public class LineageCountStatistic extends Statistic.Abstract{

	public static final String LINEAGE_COUNT_STATISTIC = "lineageCountStatistic";
	
	private OldAbstractCoalescentLikelihood acl;
	
	public LineageCountStatistic(OldAbstractCoalescentLikelihood acl){
		this.acl = acl;
	}
	
	public int getDimension() {
		return acl.getIntervalCount();
	}

	public double getStatisticValue(int dim) {
			return acl.getLineageCount(dim);
	}

	public static XMLObjectParser PARSER = new AbstractXMLObjectParser(){

		public String getParserDescription() {
			return "";
		}

		public Class getReturnType() {
			return LineageCountStatistic.class;
		}

		public XMLSyntaxRule[] getSyntaxRules() {
			return new XMLSyntaxRule[]{
					new ElementRule(OldAbstractCoalescentLikelihood.class),
			};
		}

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {
			return new LineageCountStatistic(
					(OldAbstractCoalescentLikelihood)xo.getChild(OldAbstractCoalescentLikelihood.class));
		}

		public String getParserName() {
			return LINEAGE_COUNT_STATISTIC;
		}
		
	};	
}

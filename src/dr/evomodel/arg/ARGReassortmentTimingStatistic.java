package dr.evomodel.arg;

import java.util.ArrayList;
import java.util.Collections;
import java.util.logging.Logger;

import dr.evomodel.arg.ARGModel.Node;

import dr.inference.model.Statistic;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLObjectParser;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

public class ARGReassortmentTimingStatistic extends Statistic.Abstract{

	private int dimension;
	private ARGModel arg;
	
	public static final String ARG_TIMING_STATISTIC = "argTimingStatistic";
	public static final String DIMENSION = "dimension";  //TODO This is probably somewhere else in BEAST.
	
	public ARGReassortmentTimingStatistic(String name, int dim, ARGModel arg){
		super(name);
		
		this.dimension = dim;
		this.arg = arg;
	}
	
	public int getDimension() {
		return dimension;
	}

	public double getStatisticValue(int dim) {
		if(arg.getReassortmentNodeCount() == dimension){
			ArrayList<Double> reassortmentHeights = new ArrayList<Double>();
			
			for(int i = 0; i < arg.getInternalNodeCount(); i++){
				Node x = (Node)arg.getNode(i);
				if(x.isReassortment()){
					reassortmentHeights.add(x.getHeight());
				}
			}
			
			Collections.sort(reassortmentHeights);
			return reassortmentHeights.get(dim);
		}
		return Double.NaN;
	}
	
	public static XMLObjectParser PARSER = new AbstractXMLObjectParser(){

		public String getParserDescription() {
			// TODO Auto-generated method stub
			return "";
		}

		public Class getReturnType() {
			// TODO Auto-generated method stub
			return ARGReassortmentTimingStatistic.class;
		}

		@Override
		public XMLSyntaxRule[] getSyntaxRules() {
			return new XMLSyntaxRule[]{
				new ElementRule(ARGModel.class,false),
				AttributeRule.newIntegerRule(DIMENSION,false),
				AttributeRule.newStringRule(NAME,true),
			};
		}

		@Override
		public Object parseXMLObject(XMLObject xo) throws XMLParseException {
			String name = xo.getAttribute(NAME, "");
			int dim = xo.getIntegerAttribute(DIMENSION);
			ARGModel arg = (ARGModel)xo.getChild(ARGModel.class);
			
			Logger.getLogger("dr.evomodel").info("Creating timing statistic of dimension " + dim);
			
			return new ARGReassortmentTimingStatistic(name,dim,arg);
		}

		public String getParserName() {
			return ARG_TIMING_STATISTIC;
		}
		
	};

	
	
}


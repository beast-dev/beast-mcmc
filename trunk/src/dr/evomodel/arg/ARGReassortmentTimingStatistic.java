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
	public static final String NUMBER_OF_REASSORTMENTS = "reassortments";  //TODO This is probably somewhere else in BEAST.
	
	public ARGReassortmentTimingStatistic(String name, ARGModel arg){
		super(name);
		
		
		this.dimension = arg.getExternalNodeCount() + 1;
		this.arg = arg;
	}
	
	public int getDimension() {
		return dimension;
	}

	public double getStatisticValue(int dim) {
		if(arg.getReassortmentNodeCount() != 1){
			return Double.NaN;
		}
		
		if(dim == dimension - 1){
			for(int i = 0; i < arg.getNodeCount(); i++){
				Node x = (Node)arg.getNode(i);
						
				if(x.isReassortment()){
					return (x.getHeight());
				}
				
				
			}
		}
		
		ArrayList<Double> reassortmentHeights = new ArrayList<Double>();
		
		for(int i = 0; i < arg.getInternalNodeCount(); i++){
			
			reassortmentHeights.add( ((Node)arg.getInternalNode(i)).getHeight() );
		}
		
		Collections.sort(reassortmentHeights);
		
		
		return reassortmentHeights.get(dim);
	}
	
	public static XMLObjectParser PARSER = new AbstractXMLObjectParser(){

		public String getParserDescription() {
			
			return "";
		}

		public Class getReturnType() {
			
			return ARGReassortmentTimingStatistic.class;
		}

		
		public XMLSyntaxRule[] getSyntaxRules() {
			return new XMLSyntaxRule[]{
				new ElementRule(ARGModel.class,false),
				AttributeRule.newStringRule(NAME,true),
			};
		}

		
		public Object parseXMLObject(XMLObject xo) throws XMLParseException {
			String name = xo.getId();
//			int dim = xo.getIntegerAttribute(DIMENSION);
			ARGModel arg = (ARGModel)xo.getChild(ARGModel.class);
			
			Logger.getLogger("dr.evomodel").info("Creating timing statistic");
			
			return new ARGReassortmentTimingStatistic(name,arg);
		}

		public String getParserName() {
			return ARG_TIMING_STATISTIC;
		}
		
	};

	
	
}


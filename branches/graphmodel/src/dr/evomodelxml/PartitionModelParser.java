package dr.evomodelxml;

import java.util.ArrayList;
import java.util.List;

import dr.evolution.alignment.SiteList;
import dr.evomodel.graph.PartitionModel;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

public class PartitionModelParser extends AbstractXMLObjectParser{

	public PartitionModelParser(){
		rules = new XMLSyntaxRule[]{
			new ElementRule(SiteList.class,1,Integer.MAX_VALUE),	
		};
	}
	
	public String getParserDescription() {
		return "A partition model object";
	}

	public Class getReturnType() {
		return PartitionModel.class;
	}

	public XMLSyntaxRule[] getSyntaxRules() {
		return rules;
	}

	public Object parseXMLObject(XMLObject xo) throws XMLParseException {
		List<SiteList> listOfSiteLists = new ArrayList<SiteList>();
		
		for(int i = 0; i < xo.getChildCount(); i++){
			XMLObject cxo = (XMLObject) xo.getChild(i);
			
			listOfSiteLists.add( (SiteList) cxo);
		}
				
		PartitionModel partitionModel = new PartitionModel(listOfSiteLists);
		
		return partitionModel;
	}

	public String getParserName() {
		return PartitionModel.PARTITION_MODEL;
	}

	private final XMLSyntaxRule[] rules;
	
}

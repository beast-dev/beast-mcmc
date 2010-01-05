package dr.evomodelxml;

import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

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
		
		Logger.getLogger("dr.evomodel").info("Creating a partition model, '" + xo.getId() + "', using alignments");
		for(int i = 0; i < xo.getChildCount(); i++){
			SiteList cxo = (SiteList) xo.getChild(i);
			
			listOfSiteLists.add(cxo);

			Logger.getLogger("dr.evomodel").info("\t" + cxo.getId() + "");
		}
		
		PartitionModel partitionModel = new PartitionModel(listOfSiteLists);
		
		return partitionModel;
	}

	public String getParserName() {
		return PartitionModel.PARTITION_MODEL;
	}

	private final XMLSyntaxRule[] rules;
	
}

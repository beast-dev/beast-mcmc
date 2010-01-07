package dr.evomodelxml;

import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

import dr.evolution.alignment.SiteList;
import dr.evomodel.graph.PartitionModel;
import dr.inference.model.Model;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

public class PartitionModelParser extends AbstractXMLObjectParser{

	public static final String PARTITION = "partition";
	
	public PartitionModelParser(){
		rules = new XMLSyntaxRule[]{
			
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
			XMLObject cxo = (XMLObject) xo.getChild(i);
			
			SiteList siteList = (SiteList) cxo.getChild(SiteList.class);
			
			listOfSiteLists.add(siteList);

			Logger.getLogger("dr.evomodel").info("\t" + siteList.getId() + "");
		}
		
		PartitionModel partitionModel = new PartitionModel(listOfSiteLists);
		
		for(int i = 0; i < xo.getChildCount(); i++){
			XMLObject cxo = (XMLObject) xo.getChild(i);
			
			for(int j = 0; j < cxo.getChildCount(); j++){
				if(cxo.getChild(j) instanceof Model){
					partitionModel.addModelToPartition(partitionModel.getPartition(i), 
							(Model) cxo.getChild(j));
				}else if(cxo.getChild(j) instanceof SiteList){
					//don't do anything
				}else{
					throw new XMLParseException("Something is wrong");
				}
			}
		}
		
		return partitionModel;
	}

	public String getParserName() {
		return PartitionModel.PARTITION_MODEL;
	}

	private final XMLSyntaxRule[] rules;
	
}

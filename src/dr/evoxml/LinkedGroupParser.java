package dr.evoxml;

import dr.evolution.LinkedGroup;
import dr.evolution.util.Taxa;
import dr.evolution.util.TaxonList;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

/**
 * @author Aaron Darling
 */
public class LinkedGroupParser extends AbstractXMLObjectParser {

	@Override
	public String getParserDescription() {
		return "A group of metagenome reads linked with some probability";
	}

	@Override
	public Class getReturnType() {
		return LinkedGroup.class;
	}

	@Override
	public XMLSyntaxRule[] getSyntaxRules() {
		return rules;
	}

	@Override
	public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        TaxonList taxa = null;
        double linkageProbability = 0.9999999999;
        if (xo.hasAttribute("probability")) {
            linkageProbability = xo.getDoubleAttribute("probability");
        }
    	taxa = (TaxonList)xo.getChild(TaxonList.class);
    	LinkedGroup lg = new LinkedGroup(taxa, linkageProbability);
		return lg;
	}

	@Override
	public String getParserName() {
		return "LinkedGroup";
	}
	
	private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
			new ElementRule(Taxa.class),
	        AttributeRule.newDoubleRule("probability", true, "the probability that the group of reads are linked to each other"),
    };

}

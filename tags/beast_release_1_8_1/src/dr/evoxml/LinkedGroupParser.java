package dr.evoxml;

import dr.evolution.LinkedGroup;
import dr.evolution.util.Taxa;
import dr.evolution.util.TaxonList;
import dr.xml.*;

/**
 * @author Aaron Darling
 */
public class LinkedGroupParser extends AbstractXMLObjectParser {


	public String getParserDescription() {
		return "A group of metagenome reads linked with some probability";
	}

	
	public Class getReturnType() {
		return LinkedGroup.class;
	}


	public XMLSyntaxRule[] getSyntaxRules() {
		return rules;
	}


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

	public String getParserName() {
		return "LinkedGroup";
	}
	
	private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
			new ElementRule(Taxa.class),
	        AttributeRule.newDoubleRule("probability", true, "the probability that the group of reads are linked to each other"),
    };

}

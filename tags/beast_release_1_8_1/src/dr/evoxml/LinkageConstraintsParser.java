package dr.evoxml;

import dr.evolution.LinkageConstraints;
import dr.evolution.LinkedGroup;
import dr.xml.*;

import java.util.ArrayList;

/**
 * @author Aaron Darling (koadman)
 */
public class LinkageConstraintsParser extends AbstractXMLObjectParser {

	
	public String getParserDescription() {
		return "Data representing metagenome reads that are linked by mate-pair, strobe, or other information";
	}


	public Class getReturnType() {
		return LinkageConstraints.class;
	}


	public XMLSyntaxRule[] getSyntaxRules() {
		return rules;
	}


	public Object parseXMLObject(XMLObject xo) throws XMLParseException {
    	ArrayList<LinkedGroup> groups = new ArrayList<LinkedGroup>();
        for (int i = 0; i < xo.getChildCount(); i++) {
            Object child = xo.getChild(i);
            if (child instanceof LinkedGroup) {
            	groups.add((LinkedGroup)child);
            }
    	}
    	LinkageConstraints lc = new LinkageConstraints(groups);
		return lc;
	}


	public String getParserName() {
		return "LinkageConstraints";
	}

	private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
			new ElementRule(LinkedGroup.class, 1, Integer.MAX_VALUE),
    };

}

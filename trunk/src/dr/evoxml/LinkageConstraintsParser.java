package dr.evoxml;

import java.util.ArrayList;

import dr.evolution.LinkageConstraints;
import dr.evolution.LinkedGroup;
import dr.evolution.alignment.Alignment;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;
import dr.xml.XORRule;

/**
 * @author Aaron Darling (koadman)
 */
public class LinkageConstraintsParser extends AbstractXMLObjectParser {

	@Override
	public String getParserDescription() {
		return "Data representing metagenome reads that are linked by mate-pair, strobe, or other information";
	}

	@Override
	public Class getReturnType() {
		return LinkageConstraints.class;
	}

	@Override
	public XMLSyntaxRule[] getSyntaxRules() {
		return rules;
	}

	@Override
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

	@Override
	public String getParserName() {
		return "LinkageConstraints";
	}

	private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
			new ElementRule(LinkedGroup.class, 1, Integer.MAX_VALUE),
    };

}

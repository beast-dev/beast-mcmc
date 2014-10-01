package dr.evoxml;

import dr.evolution.LinkageConstraints;
import dr.evolution.MetagenomeData;
import dr.evolution.alignment.Alignment;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxa;
import dr.evolution.util.TaxonList;
import dr.xml.*;

/**
 * @author Aaron Darling
 */
public class MetagenomeDataParser extends AbstractXMLObjectParser {

	@Override
	public String getParserDescription() {
		return "Data representing metagenome reads aligned to reference sequences and a reference tree";
	}

	@Override
	public Class getReturnType() {
		return MetagenomeData.class;
	}

	@Override
	public XMLSyntaxRule[] getSyntaxRules() {
		return rules;
	}

	@Override
	public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        TaxonList taxa = null;
    	taxa = (TaxonList)xo.getChild(TaxonList.class);
        if(taxa==null)
        	taxa = (Tree)xo.getChild(Tree.class);

        Alignment alignment = (Alignment) xo.getChild(Alignment.class);
        LinkageConstraints lc = (LinkageConstraints)xo.getChild(LinkageConstraints.class);
        
        boolean fixedReferenceTree = false;
        if (xo.hasAttribute("fixedReferenceTree")) {
        	fixedReferenceTree = xo.getBooleanAttribute("fixedReferenceTree");
        }

        MetagenomeData md = new MetagenomeData(taxa, alignment, lc, fixedReferenceTree); 
        return md;
	}

	
	public String getParserName() {
		return "MetagenomeData";
	}

	private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
			new XORRule(new ElementRule(Taxa.class),
						new ElementRule(Tree.class)),
            new ElementRule(Alignment.class),
            new ElementRule(LinkageConstraints.class, true),	// optional element
	        AttributeRule.newBooleanRule("fixedReferenceTree", true, "Whether the reference tree should be of fixed topology"),
    };

}

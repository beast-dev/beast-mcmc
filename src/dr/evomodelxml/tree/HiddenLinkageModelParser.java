package dr.evomodelxml.tree;

import dr.evolution.MetagenomeData;
import dr.evomodel.tree.HiddenLinkageModel;
import dr.evomodel.tree.TreeModel;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.ElementRule;
import dr.xml.StringAttributeRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

/**
 * @author Aaron Darling (koadman)
 */
public class HiddenLinkageModelParser extends AbstractXMLObjectParser {

    public static final String LINKAGE_GROUP_COUNT = "linkageGroupCount";

    public static final String NAME = "HiddenLinkageModel";

    @Override
	public String getParserDescription() {
		return "A model to describe missing information about linkage among several reads from a metagenome";
	}

	@Override
	public Class getReturnType() {
		return HiddenLinkageModel.class;
	}

	@Override
	public XMLSyntaxRule[] getSyntaxRules() {
		return rules;
	}

	@Override
	public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        String linkageGroupCount = xo.getAttribute(LINKAGE_GROUP_COUNT, xo.getId());
        MetagenomeData data = (MetagenomeData)xo.getChild(MetagenomeData.class);

        int tc = Integer.parseInt(linkageGroupCount);
        return new HiddenLinkageModel(tc, data);
	}

	@Override
	public String getParserName() {
		return NAME;
	}

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new StringAttributeRule(LINKAGE_GROUP_COUNT, "The number of hidden lineages", true),
            new ElementRule(MetagenomeData.class)
    };
}

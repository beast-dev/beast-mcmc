package dr.evomodelxml.tree;

import dr.evolution.MetagenomeData;
import dr.evomodel.tree.HiddenLinkageModel;
import dr.xml.*;

/**
 * @author Aaron Darling (koadman)
 */
public class HiddenLinkageModelParser extends AbstractXMLObjectParser {

    public static final String LINKAGE_GROUP_COUNT = "linkageGroupCount";

    public static final String NAME = "HiddenLinkageModel";

    
	public String getParserDescription() {
		return "A model to describe missing information about linkage among several reads from a metagenome";
	}


	public Class getReturnType() {
		return HiddenLinkageModel.class;
	}


	public XMLSyntaxRule[] getSyntaxRules() {
		return rules;
	}


	public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        String linkageGroupCount = xo.getAttribute(LINKAGE_GROUP_COUNT, xo.getId());
        MetagenomeData data = (MetagenomeData)xo.getChild(MetagenomeData.class);

        int tc = Integer.parseInt(linkageGroupCount);
        return new HiddenLinkageModel(tc, data);
	}


	public String getParserName() {
		return NAME;
	}

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new StringAttributeRule(LINKAGE_GROUP_COUNT, "The number of hidden lineages", true),
            new ElementRule(MetagenomeData.class)
    };
}

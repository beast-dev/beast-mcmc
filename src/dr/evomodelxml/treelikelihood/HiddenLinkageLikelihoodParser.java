package dr.evomodelxml.treelikelihood;

import dr.evomodel.tree.HiddenLinkageModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treelikelihood.HiddenLinkageLikelihood;
import dr.evomodel.treelikelihood.TreeLikelihood;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

/**
 * @author Aaron Darling (koadman)
 */
public class HiddenLinkageLikelihoodParser extends AbstractXMLObjectParser {

	@Override
	public String getParserDescription() {
		return "A likelihood calculator for hidden linkage among metagenomic reads";
	}

	@Override
	public Class getReturnType() {
		return HiddenLinkageLikelihood.class;
	}

	@Override
	public XMLSyntaxRule[] getSyntaxRules() {
		return rules;
	}

	@Override
	public Object parseXMLObject(XMLObject xo) throws XMLParseException {
		HiddenLinkageModel hlm = (HiddenLinkageModel) xo.getChild(HiddenLinkageModel.class);
		TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);
        return new HiddenLinkageLikelihood(hlm, tree);
	}

	@Override
	public String getParserName() {
		return "HiddenLinkageLikelihood";
	}

	private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(HiddenLinkageModel.class),
            new ElementRule(TreeModel.class),
    };

}

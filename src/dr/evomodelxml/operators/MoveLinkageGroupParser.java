package dr.evomodelxml.operators;

import dr.evolution.MetagenomeData;
import dr.evomodel.operators.MoveLinkageGroup;
import dr.evomodel.operators.WilsonBalding;
import dr.evomodel.tree.HiddenLinkageModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.operators.MCMCOperator;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.ElementRule;
import dr.xml.StringAttributeRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

/**
 * @author Aaron Darling (koadman)
 */
public class MoveLinkageGroupParser extends AbstractXMLObjectParser {
    public static final String MOVE_LINKAGE_GROUP = "moveLinkageGroup";

	public String getParserDescription() {
		return "Operator to reassign metagenomic reads from one linkage group to another";
	}

	public Class getReturnType() {
		return MoveLinkageGroup.class;
	}

	public XMLSyntaxRule[] getSyntaxRules() {
		return rules;
	}

	public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        final double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
        final HiddenLinkageModel hlm = (HiddenLinkageModel) xo.getChild(HiddenLinkageModel.class);
        return new MoveLinkageGroup(hlm, weight);
	}

	public String getParserName() {
		return MOVE_LINKAGE_GROUP;
	}

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new StringAttributeRule(MCMCOperator.WEIGHT, "Weight of the move", true),
            new ElementRule(HiddenLinkageModel.class)
    };
}

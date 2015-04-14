package dr.app.beagle.evomodel.branchmodel.lineagespecific;

import dr.app.beagle.evomodel.branchmodel.BranchModel;
import dr.evomodel.tree.TreeModel;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

/**
 * @author Filip Bielejec
 * @version $Id$
 * 
 */
public class BranchSpecificTraitParser extends AbstractXMLObjectParser {

	public static final String BRANCH_SPECIFIC_TRAIT = "branchSpecificTrait";

	@Override
	public String getParserName() {
		return BRANCH_SPECIFIC_TRAIT;
	}

	@Override
	public Object parseXMLObject(XMLObject xo) throws XMLParseException {

		BranchModel branchModel = (BranchModel) xo.getChild(BranchModel.class);
//		CompoundParameter parameter = (CompoundParameter) xo.getChild(CompoundParameter.class);
        TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
		
		return new BranchSpecificTrait(treeModel, branchModel, xo.getId());
	}

	@Override
	public XMLSyntaxRule[] getSyntaxRules() {
		return new XMLSyntaxRule[] {

		new ElementRule(BranchModel.class, false), //
				new ElementRule(TreeModel.class, false), //

		};
	}

	@Override
	public String getParserDescription() {
		return BRANCH_SPECIFIC_TRAIT;
	}

	@Override
	public Class getReturnType() {
		return BranchSpecificTrait.class;
	}

}// END: class

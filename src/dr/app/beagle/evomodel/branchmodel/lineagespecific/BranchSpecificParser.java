package dr.app.beagle.evomodel.branchmodel.lineagespecific;

import dr.app.beagle.evomodel.substmodel.SubstitutionModel;
import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.CountableBranchCategoryProvider;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.tree.TreeModelParser;
import dr.inference.model.Parameter;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

public class BranchSpecificParser extends AbstractXMLObjectParser {

	public static final String BRANCH_SPECIFIC = "branchSpecificBranchRates";
	public static final String RATES = "rates";
	
	@Override
	public String getParserName() {
		return BRANCH_SPECIFIC;
	}

	@Override
	public Object parseXMLObject(XMLObject xo) throws XMLParseException {

		TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
		
		CountableBranchCategoryProvider rateCategories = null;
        Parameter ratesParameter = null;
		
		 if (xo.hasChildNamed(RATES )) {
			 
			 ratesParameter = (Parameter) xo.getElementFirstChild(RATES);
			 
			 rateCategories = new CountableBranchCategoryProvider.IndependentBranchCategoryModel(treeModel, ratesParameter);
			 
		 }
		
			SubstitutionModel codonModel = (SubstitutionModel) xo.getChild(SubstitutionModel.class);
		 
		return new BranchSpecific(treeModel, codonModel, rateCategories);
	}

	@Override
	public XMLSyntaxRule[] getSyntaxRules() {
		return null;
	}

	@Override
	public String getParserDescription() {
		return BRANCH_SPECIFIC;
	}

	@Override
	public Class getReturnType() {
		return BranchSpecific.class;
	}

}//END: class

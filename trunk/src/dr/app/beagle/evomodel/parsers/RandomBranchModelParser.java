package dr.app.beagle.evomodel.parsers;

import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

import dr.app.beagle.evomodel.branchmodel.RandomBranchModel;
import dr.app.beagle.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.tree.TreeModel;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

public class RandomBranchModelParser extends AbstractXMLObjectParser {

	 public static final String MODELS = "models";
	
	@Override
	public String getParserName() {
		return RandomBranchModel.RANDOM_BRANCH_MODEL;
	}

	@Override
	public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Logger.getLogger("dr.evomodel").info("Using random assignment branch model.");
        TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);

	      XMLObject cxo = xo.getChild(MODELS);
	      List<SubstitutionModel> substitutionModels = new ArrayList<SubstitutionModel>();
	      for (int i = 0; i < cxo.getChildCount(); i++) {

				SubstitutionModel substModel = (SubstitutionModel) cxo.getChild(i);
				substitutionModels.add(substModel);
				
			}//END: models loop
		
		return new RandomBranchModel(treeModel, substitutionModels);
	}//END: parseXMLObject

	@Override
	public XMLSyntaxRule[] getSyntaxRules() {
		
		return new XMLSyntaxRule[]{
				
				new ElementRule(TreeModel.class, false), //
                new ElementRule(MODELS,
                        new XMLSyntaxRule[] {
                                new ElementRule(SubstitutionModel.class, 1, Integer.MAX_VALUE),
                        }
                )             
                
		};
	}//END: XMLSyntaxRule

	@Override
	public String getParserDescription() {
		return "This element provides a branch model which randomly assigns " +
				"substitution models to branches on the tree by sampling " +
				"with replacement from the provided list of substitution models. ";
	}

	@Override
	public Class getReturnType() {
		return RandomBranchModel.class;
	}


}//END: class

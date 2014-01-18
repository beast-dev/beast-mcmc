package dr.app.beagle.evomodel.branchmodel.lineagespecific;

import java.util.ArrayList;
import java.util.List;

import dr.app.beagle.evomodel.substmodel.FrequencyModel;
import dr.app.beagle.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.branchratemodel.CountableBranchCategoryProvider;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.AbstractModel;
import dr.inference.model.Parameter;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

public class BranchSpecificParser extends AbstractXMLObjectParser {

	 public static final String BRANCH_SPECIFIC = "branchSpecific";
	 public static final String MODELS = "models";
	 public static final String CATEGORIES = "categories";
	
	@Override
	public String getParserName() {
		return BRANCH_SPECIFIC;
	}

	@Override
	public Object parseXMLObject(XMLObject xo) throws XMLParseException {

		TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
		FrequencyModel rootFrequencyModel = (FrequencyModel) xo.getChild(FrequencyModel.class);
		
	      List<SubstitutionModel> substitutionModels = new ArrayList<SubstitutionModel>();
	      
	      XMLObject cxo = xo.getChild(MODELS);
			for (int i = 0; i < cxo.getChildCount(); i++) {

				SubstitutionModel substModel = (SubstitutionModel) cxo.getChild(i);
				substitutionModels.add(substModel);
				
			}//END: models loop
		
			// TODO: check if categories numbering starts from zero
			Parameter categories = (Parameter) xo.getElementFirstChild(CATEGORIES); 
			
			CountableBranchCategoryProvider.CladeBranchCategoryModel provider = new CountableBranchCategoryProvider.CladeBranchCategoryModel(treeModel, categories, false);
			
		return new BranchSpecific(treeModel, rootFrequencyModel, substitutionModels, provider, categories);
	}

	@Override
	public XMLSyntaxRule[] getSyntaxRules() {
		return new XMLSyntaxRule[]{
				
				new ElementRule(TreeModel.class, false), //
				new ElementRule(FrequencyModel.class, false), //
                new ElementRule(MODELS,
                        new XMLSyntaxRule[] {
                                new ElementRule(AbstractModel.class, 1, Integer.MAX_VALUE),
                        }
                ), //
                new ElementRule(CATEGORIES,
                        new XMLSyntaxRule[] {
                                new ElementRule(Parameter.class, 1, 1),
                        }
                ) //              
                
                
		};
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

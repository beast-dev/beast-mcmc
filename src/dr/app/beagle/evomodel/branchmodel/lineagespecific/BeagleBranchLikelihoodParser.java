package dr.app.beagle.evomodel.branchmodel.lineagespecific;

import java.util.ArrayList;
import java.util.List;

import dr.app.beagle.evomodel.branchmodel.lineagespecific.BeagleBranchLikelihood;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

public class BeagleBranchLikelihoodParser extends AbstractXMLObjectParser {

	public static final String BEAGLE_BRANCH_LIKELIHOODS = "beagleBranchLikelihood";
	
	public static final String TREE_LIKELIHOODS = "treelLikelihoods";
	
	@Override
	public String getParserName() {
		return BEAGLE_BRANCH_LIKELIHOODS;
	}

	@Override
	public Object parseXMLObject(XMLObject xo) throws XMLParseException {

	       TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
		
	       Parameter zParameter = (Parameter) xo.getElementFirstChild(  DirichletProcessPriorParser.CATEGORIES);
	       
	        List<Likelihood> treeLikelihoods = new ArrayList<Likelihood>();
	       
	        XMLObject cxo = (XMLObject) xo.getChild(TREE_LIKELIHOODS);
	        for (int i = 0; i < cxo.getChildCount(); i++) {
	        
	        	Likelihood likelihood = (Likelihood) cxo.getChild(i);
	        	treeLikelihoods.add(likelihood);
	        }
	        
	        
		return new BeagleBranchLikelihood(treeModel, treeLikelihoods, zParameter);
	}

	@Override
	public XMLSyntaxRule[] getSyntaxRules() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public String getParserDescription() {
		return BEAGLE_BRANCH_LIKELIHOODS;
	}

	@Override
	public Class getReturnType() {
		return BeagleBranchLikelihood.class;
	}


}

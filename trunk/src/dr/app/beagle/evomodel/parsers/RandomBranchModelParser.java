package dr.app.beagle.evomodel.parsers;

import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

import dr.app.beagle.evomodel.branchmodel.RandomBranchModel;
import dr.app.beagle.evomodel.substmodel.GY94CodonModel;
import dr.app.beagle.evomodel.substmodel.SubstitutionModel;
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

public class RandomBranchModelParser extends AbstractXMLObjectParser {

	 public static final String BASE_MODEL = "baseSubstitutionModel";
	
	@Override
	public String getParserName() {
		return RandomBranchModel.RANDOM_BRANCH_MODEL;
	}

	@Override
	public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Logger.getLogger("dr.evomodel").info("Using random assignment branch model.");
        
        TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
        GY94CodonModel baseSubstitutionModel = (GY94CodonModel) xo .getElementFirstChild(BASE_MODEL);
        
		return new RandomBranchModel(treeModel, baseSubstitutionModel);
	}//END: parseXMLObject

	@Override
	public XMLSyntaxRule[] getSyntaxRules() {

		return new XMLSyntaxRule[] {

		new ElementRule(TreeModel.class, false), //
        new ElementRule(BASE_MODEL,
                new XMLSyntaxRule[]{new ElementRule(SubstitutionModel.class, 1, 1)})

		};
	}// END: XMLSyntaxRule

	@Override
	public String getParserDescription() {
		return RandomBranchModel.RANDOM_BRANCH_MODEL;
	}

	@Override
	public Class getReturnType() {
		return RandomBranchModel.class;
	}


}//END: class

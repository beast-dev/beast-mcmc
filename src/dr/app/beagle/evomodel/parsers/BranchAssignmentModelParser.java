package dr.app.beagle.evomodel.parsers;

import java.util.LinkedHashMap;
import java.util.logging.Logger;

import dr.app.beagle.evomodel.branchmodel.BranchAssignmentModel;
import dr.app.beagle.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.tree.TreeModel;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

public class BranchAssignmentModelParser extends AbstractXMLObjectParser {

	public static final String ANNOTATION = "annotation";
	public static final String ANNOTATION_VALUE = "annotationValue";
	public static final String BASE_MODEL = "baseModel";
	public static final String ASSIGNMENT = "assignment";

	@Override
	public String getParserName() {
		return BranchAssignmentModel.BRANCH_ASSIGNMENT_MODEL;
	}

	@Override
	public Object parseXMLObject(XMLObject xo) throws XMLParseException {

		Logger.getLogger("dr.evomodel").info(
				"Using branch assignment branch model.");
		TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
		String annotation = xo.getStringAttribute(ANNOTATION);

		LinkedHashMap<Integer, SubstitutionModel> modelIndexMap = new LinkedHashMap<Integer, SubstitutionModel>();

		for (int i = 0; i < xo.getChildCount(); i++) {
			if (xo.getChild(i) instanceof XMLObject) {
				XMLObject xoc = (XMLObject) xo.getChild(i);
				if (xoc.getName().equals(ASSIGNMENT)) {

					Integer index = null;
					if (xoc.hasAttribute(ANNOTATION_VALUE)) {
						index = xoc.getIntegerAttribute(ANNOTATION_VALUE);
					}

					SubstitutionModel model = (SubstitutionModel) xoc
							.getChild(SubstitutionModel.class);

					modelIndexMap.put(index, model);

				}
			}
		}

		SubstitutionModel baseModel = (SubstitutionModel) xo
				.getElementFirstChild(BASE_MODEL);

		return new BranchAssignmentModel(treeModel, annotation, modelIndexMap,
				baseModel);
	}// END: parseXMLObject

	@Override
	public XMLSyntaxRule[] getSyntaxRules() {

		return new XMLSyntaxRule[] {

		new ElementRule(TreeModel.class, false), //
				new ElementRule(ASSIGNMENT, //
						new XMLSyntaxRule[] { //
								AttributeRule.newIntegerRule(ANNOTATION_VALUE,
										false), //
								new ElementRule(SubstitutionModel.class, false), // model
						},//
						1, Integer.MAX_VALUE), //
				new ElementRule(BASE_MODEL, //
						new XMLSyntaxRule[] { //
						new ElementRule(SubstitutionModel.class, false) //
						}, false) // base model

		};
	}// END: getSyntaxRules

	@Override
	public String getParserDescription() {
		return "This element provides a branch model which assigns "
				+ "substitution models to branches on the tree ";
	}

	@SuppressWarnings("rawtypes")
	@Override
	public Class getReturnType() {
		return BranchAssignmentModelParser.class;
	}

}// END: class

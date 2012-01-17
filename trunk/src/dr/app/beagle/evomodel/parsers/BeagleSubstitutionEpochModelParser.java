package dr.app.beagle.evomodel.parsers;

import java.util.ArrayList;
import java.util.List;

import dr.app.beagle.evomodel.sitemodel.EpochBranchSubstitutionModel;
import dr.app.beagle.evomodel.substmodel.FrequencyModel;
import dr.app.beagle.evomodel.substmodel.SubstitutionModel;
import dr.evolution.datatype.DataType;
import dr.inference.model.Parameter;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

public class BeagleSubstitutionEpochModelParser extends AbstractXMLObjectParser {

	public static final String SUBSTITUTION_EPOCH_MODEL = "beagleSubstitutionEpochModel";
	public static final String MODELS = "models";
	public static final String TRANSITION_TIMES = "transitionTimes";

	@Override
	public Object parseXMLObject(XMLObject xo) throws XMLParseException {

		DataType dataType = null;
		List<FrequencyModel> frequencyModelList = new ArrayList<FrequencyModel>();
		List<SubstitutionModel> substModelList = new ArrayList<SubstitutionModel>();
		XMLObject cxo = xo.getChild(MODELS);
		
		for (int i = 0; i < cxo.getChildCount(); i++) {
			SubstitutionModel model = (SubstitutionModel) cxo.getChild(i);

			if (dataType == null) {
				dataType = model.getDataType();
			} else if (dataType != model.getDataType())
				throw new XMLParseException(
						"Substitution models across epoches must use the same data type.");

			if (frequencyModelList.size() == 0) {

				frequencyModelList.add(model.getFrequencyModel()); // model.getFrequencyModel();

			} else if (frequencyModelList.get(0) != model.getFrequencyModel())
				throw new XMLParseException(
						"Substitution models across epoches must currently use the same frequency model.\n Harass Marc to fix this.");

			substModelList.add(model);
		}

		Parameter epochTransitionTimes = (Parameter) xo
				.getChild(Parameter.class);

		if (epochTransitionTimes.getDimension() != substModelList.size() - 1) {
			throw new XMLParseException(
					"# of transition times must equal # of substitution models - 1\n"
							+ epochTransitionTimes.getDimension() + "\n"
							+ substModelList.size());
		}

		
		//TODO: sort in increasing order
		
		
		return new EpochBranchSubstitutionModel(substModelList,
				frequencyModelList, epochTransitionTimes);
	}// END: parseXMLObject

	@Override
	public XMLSyntaxRule[] getSyntaxRules() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public String getParserDescription() {
		return null;
	}

	@Override
	public Class getReturnType() {
		return EpochBranchSubstitutionModel.class;
	}

	@Override
	public String getParserName() {
		return SUBSTITUTION_EPOCH_MODEL;
	}

}// END: class

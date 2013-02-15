package dr.evomodelxml.substmodel;

import dr.evolution.datatype.DataType;
import dr.evomodel.substmodel.AbstractSubstitutionModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.SubstitutionEpochModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

/**
 */
public class SubstitutionEpochModelParser extends AbstractXMLObjectParser {

    public static final String SUBSTITUTION_EPOCH_MODEL = "substitutionEpochModel";
    public static final String MODELS = "models";
    public static final String TRANSITION_TIMES = "transitionTimes";

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        DataType dataType = null;
        FrequencyModel freqModel = null;
        List<SubstitutionModel> modelList = new ArrayList<SubstitutionModel>();
        XMLObject cxo = xo.getChild(MODELS);
        for (int i = 0; i < cxo.getChildCount(); i++) {
            SubstitutionModel model = (SubstitutionModel) cxo.getChild(i);

            if (dataType == null) {
                dataType = model.getDataType();
            } else if (dataType != model.getDataType())
                throw new XMLParseException("Substitution models across epoches must use the same data type.");

            if (freqModel == null) {
                freqModel = model.getFrequencyModel();
            } else if (freqModel != model.getFrequencyModel())
                throw new XMLParseException("Substitution models across epoches must currently use the same frequency model.\nHarass Marc to fix this.");

            modelList.add(model);
        }

        Parameter transitionTimes = (Parameter) xo.getChild(Parameter.class);

        if (transitionTimes.getDimension() != modelList.size() - 1) {
            throw new XMLParseException("# of transition times must equal # of substitution models - 1\n" + transitionTimes.getDimension() + "\n" + modelList.size());
        }

        return new SubstitutionEpochModel(SUBSTITUTION_EPOCH_MODEL, modelList, transitionTimes, dataType, freqModel);
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                new ElementRule(MODELS,
                        new XMLSyntaxRule[]{
                                new ElementRule(AbstractSubstitutionModel.class, 1, Integer.MAX_VALUE),
                        }),
                new ElementRule(Parameter.class),
        };
    }

    public String getParserDescription() {
        return null;
    }

    public Class getReturnType() {
        return SubstitutionEpochModel.class;
    }

    public String getParserName() {
        return SUBSTITUTION_EPOCH_MODEL;
    }
}

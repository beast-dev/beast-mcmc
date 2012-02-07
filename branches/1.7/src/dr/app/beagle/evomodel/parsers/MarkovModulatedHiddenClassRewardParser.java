package dr.app.beagle.evomodel.parsers;

import dr.app.beagle.evomodel.substmodel.MarkovModulatedSubstitutionModel;
import dr.evolution.datatype.HiddenDataType;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * @author Marc A. Suchard
 */

public class MarkovModulatedHiddenClassRewardParser extends AbstractXMLObjectParser {

    public static final String PARSER_NAME = "hiddenClassRewardParameter";
    public static final String CLASS_NUMBER = "class";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        MarkovModulatedSubstitutionModel substitutionModel = (MarkovModulatedSubstitutionModel) xo.getChild(MarkovModulatedSubstitutionModel.class);
        HiddenDataType hiddenDataType = (HiddenDataType) substitutionModel.getDataType();
        int classNumber = xo.getIntegerAttribute(CLASS_NUMBER);
        int hiddenClassCount = hiddenDataType.getHiddenClassCount();
        if (classNumber < 1 || classNumber > hiddenClassCount) {
            throw new XMLParseException("Invalid class number in " + xo.getId());
        }
        classNumber--; // Use zero-indexed number
        int stateCount = hiddenDataType.getStateCount() / hiddenClassCount;

        // Construct reward parameter
        Parameter parameter = new Parameter.Default(stateCount * hiddenClassCount, 0.0);
        for (int i = 0; i < stateCount; ++i) {
            parameter.setParameterValue(i + classNumber * stateCount, 1.0);
        }

        return parameter;
    }

    /**
     * @return an array of syntax rules required by this element.
     *         Order is not important.
     */
    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = {
            AttributeRule.newIntegerRule(CLASS_NUMBER),
            new ElementRule(MarkovModulatedSubstitutionModel.class),
    };

    @Override
    public String getParserDescription() {
        return "Generates a reward parameter to log hidden classes in Markov-modulated substitutionProcess";
    }

    @Override
    public Class getReturnType() {
        return Parameter.class;
    }

    /**
     * @return Parser name, which is identical to name of xml element parsed by it.
     */
    public String getParserName() {
        return PARSER_NAME;
    }
}

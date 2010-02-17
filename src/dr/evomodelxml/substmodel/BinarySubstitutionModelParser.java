package dr.evomodelxml.substmodel;

import dr.evolution.datatype.DataType;
import dr.evolution.datatype.TwoStates;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.GeneralSubstitutionModel;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * @author Michael Defoin Platel
 */
public class BinarySubstitutionModelParser extends AbstractXMLObjectParser {

    public static final String BINARY_SUBSTITUTION_MODEL = "binarySubstitutionModel";

    public String getParserName() {
        return BINARY_SUBSTITUTION_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Parameter ratesParameter;

        XMLObject cxo = xo.getChild(GeneralSubstitutionModelParser.FREQUENCIES);
        FrequencyModel freqModel = (FrequencyModel) cxo.getChild(FrequencyModel.class);

        DataType dataType = freqModel.getDataType();

        if (dataType != TwoStates.INSTANCE)
            throw new XMLParseException("Frequency model must have binary (two state) data type.");

        int relativeTo = 0;

        ratesParameter = new Parameter.Default(0);

        return new GeneralSubstitutionModel(dataType, freqModel, ratesParameter, relativeTo);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A general reversible model of sequence substitution for binary data type.";
    }

    public Class getReturnType() {
        return GeneralSubstitutionModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(GeneralSubstitutionModelParser.FREQUENCIES, FrequencyModel.class),
    };

}

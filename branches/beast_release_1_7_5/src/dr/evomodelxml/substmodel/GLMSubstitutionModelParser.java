package dr.evomodelxml.substmodel;

import dr.evolution.datatype.DataType;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.GLMSubstitutionModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evoxml.util.DataTypeUtils;
import dr.inference.distribution.GeneralizedLinearModel;
import dr.inference.distribution.LogLinearModel;
import dr.xml.*;

/**
 */
public class GLMSubstitutionModelParser extends AbstractXMLObjectParser {

    public static final String GLM_SUBSTITUTION_MODEL = "glmSubstitutionModel";


    public String getParserName() {
        return GLM_SUBSTITUTION_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        DataType dataType = DataTypeUtils.getDataType(xo);

        if (dataType == null) dataType = (DataType) xo.getChild(DataType.class);

        int rateCount = (dataType.getStateCount() - 1) * dataType.getStateCount();

        LogLinearModel glm = (LogLinearModel) xo.getChild(GeneralizedLinearModel.class);

        int length = glm.getXBeta().length;

        if (length != rateCount) {
            throw new XMLParseException("Rates parameter in " + getParserName() + " element should have " + (rateCount) + " dimensions.  However GLM dimension is " + length);
        }

        XMLObject cxo = xo.getChild(ComplexSubstitutionModelParser.ROOT_FREQUENCIES);
        FrequencyModel rootFreq = (FrequencyModel) cxo.getChild(FrequencyModel.class);

        if (dataType != rootFreq.getDataType()) {
            throw new XMLParseException("Data type of " + getParserName() + " element does not match that of its rootFrequencyModel.");
        }

        return new GLMSubstitutionModel(xo.getId(), dataType, rootFreq, glm);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A general model of sequence substitution for any data type where the rates come from the generalized linear model.";
    }

    public Class getReturnType() {
        return SubstitutionModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new XORRule(
                    new StringAttributeRule(DataType.DATA_TYPE, "The type of sequence data",
                            DataType.getRegisteredDataTypeNames(), false),
                    new ElementRule(DataType.class)
            ),
            new ElementRule(ComplexSubstitutionModelParser.ROOT_FREQUENCIES, FrequencyModel.class),
            new ElementRule(GeneralizedLinearModel.class),
    };

}

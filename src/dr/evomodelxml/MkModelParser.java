package dr.evomodelxml;

import dr.evolution.datatype.DataType;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.GeneralSubstitutionModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.logging.Logger;

/**
 * @author Alexander V Alekseyenko (alexander.alekseyenko@gmail.com)
 *         <p/>
 *         MkModelParser implements Lewis's Mk Model Syst. Biol. 50(6):913-925, 2001
 */
public class MkModelParser extends AbstractXMLObjectParser {

    public static final String MK_SUBSTITUTION_MODEL = "mkSubstitutionModel";

    public String getParserName() {
        return MK_SUBSTITUTION_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        XMLObject cxo = (XMLObject) xo.getChild(GeneralSubstitutionModel.FREQUENCIES);
        FrequencyModel freqModel = (FrequencyModel) cxo.getChild(FrequencyModel.class);

        DataType dataType = freqModel.getDataType();
        int rateCount = ((dataType.getStateCount() - 1) * dataType.getStateCount()) / 2 - 1;
        Parameter ratesParameter = new Parameter.Default(rateCount, 1.0);

        Logger.getLogger("dr.evolution").info("Creating an Mk substitution model with data type: " + dataType.getType() + "on " + dataType.getStateCount() + " states.");

        int relativeTo = 1;

        return new GeneralSubstitutionModel(dataType, freqModel, ratesParameter, relativeTo);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "An Mk model of substitution. This model can also accomodate arbitrary orderings of changes.";
    }

    public Class getReturnType() {
        return SubstitutionModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(GeneralSubstitutionModel.FREQUENCIES, FrequencyModel.class),
    };

}
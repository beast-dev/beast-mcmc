package dr.app.beagle.evomodel.parsers;

import dr.xml.*;
import dr.inference.model.Parameter;
import dr.app.beagle.evomodel.substmodel.FrequencyModel;
import dr.app.beagle.evomodel.substmodel.SubstitutionModel;
import dr.app.beagle.evomodel.substmodel.ComplexSubstitutionModel;
import dr.evolution.datatype.DataType;

import java.util.logging.Logger;

/**
 * @author Marc Suchard
 */
public class ComplexSubstitutionModelParser extends AbstractXMLObjectParser {

    public static final String COMPLEX_SUBSTITUTION_MODEL = "complexSubstitutionModel";
    public static final String DATA_TYPE = "dataType";
    public static final String RATES = "rates";
    public static final String FREQUENCIES = "frequencies";


    public String getParserName() {
        return COMPLEX_SUBSTITUTION_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Parameter ratesParameter;

        XMLObject cxo = xo.getChild(FREQUENCIES);
        FrequencyModel freqModel = (FrequencyModel) cxo.getChild(FrequencyModel.class);

        DataType dataType = freqModel.getDataType();

        cxo = xo.getChild(RATES);

        int states = dataType.getStateCount();

        Logger.getLogger("dr.app.beagle.evomodel").info("  Complex Substitution Model (stateCount=" + states + ")");

        ratesParameter = (Parameter) cxo.getChild(Parameter.class);

        int rateCount = (dataType.getStateCount() - 1) * dataType.getStateCount();

        if (ratesParameter == null) {

            if (rateCount == 1) {
                // simplest model for binary traits...
            } else {
                throw new XMLParseException("No rates parameter found in " + getParserName());
            }
        } else if (ratesParameter.getDimension() != rateCount) {
            throw new XMLParseException("Rates parameter in " + getParserName() + " element should have " + rateCount + " dimensions.");
        }

        return new ComplexSubstitutionModel(COMPLEX_SUBSTITUTION_MODEL,dataType, freqModel, ratesParameter);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A general irreversible model of sequence substitution for any data type.";
    }

    public Class getReturnType() {
        return SubstitutionModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(FREQUENCIES, FrequencyModel.class),
            new ElementRule(RATES,
                new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class, true)}
            )
    };
}

package dr.app.beagle.evomodel.parsers;

import dr.app.beagle.evomodel.substmodel.FrequencyModel;
import dr.app.beagle.evomodel.substmodel.GeneralSubstitutionModel;
import dr.app.beagle.evomodel.substmodel.SubstitutionModel;
import dr.evolution.datatype.DataType;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.logging.Logger;

/**
 * @author Marc Suchard
 * @author Andrew Rambaut
 */

public class GeneralSubstitutionModelParser extends AbstractXMLObjectParser {


    public static final String GENERAL_SUBSTITUTION_MODEL = "generalSubstitutionModel";
    public static final String DATA_TYPE = "dataType";
    public static final String RATES = "rates";
    public static final String RELATIVE_TO = "relativeTo";
    public static final String FREQUENCIES = "frequencies";


    public String getParserName() {
        return GENERAL_SUBSTITUTION_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Parameter ratesParameter;

        XMLObject cxo = xo.getChild(FREQUENCIES);
        FrequencyModel freqModel = (FrequencyModel) cxo.getChild(FrequencyModel.class);

        DataType dataType = freqModel.getDataType();

        cxo = xo.getChild(RATES);

        int states = dataType.getStateCount();

        Logger.getLogger("dr.app.beagle.evomodel").info("  General Substitution Model (stateCount=" + states + ")");

        int relativeTo = cxo.getIntegerAttribute(RELATIVE_TO) - 1;
        if (relativeTo < 0) {
            throw new XMLParseException(RELATIVE_TO + " must be 1 or greater");
        } else {
            int t = relativeTo;
            int s = states - 1;
            int row = 0;
            while (t >= s) {
                t -= s;
                s -= 1;
                row += 1;
            }
            int col = t + row + 1;

            Logger.getLogger("dr.app.beagle.evomodel").info("  Rates relative to "
                    + dataType.getCode(row) + "<->" + dataType.getCode(col));
        }

        ratesParameter = (Parameter) cxo.getChild(Parameter.class);

        int rateCount = ((dataType.getStateCount() - 1) * dataType.getStateCount()) / 2;

        if (ratesParameter == null) {

            if (rateCount == 1) {
                // simplest model for binary traits...
            } else {
                throw new XMLParseException("No rates parameter found in " + getParserName());
            }
        } else if (ratesParameter.getDimension() != rateCount - 1) {
            throw new XMLParseException("Rates parameter in " + getParserName() + " element should have " + (rateCount - 1) + " dimensions.");
        }

        return new GeneralSubstitutionModel(GENERAL_SUBSTITUTION_MODEL,dataType, freqModel, ratesParameter, relativeTo);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A general reversible model of sequence substitution for any data type.";
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
                            AttributeRule.newIntegerRule(RELATIVE_TO, false, "The index of the implicit rate (value 1.0) that all other rates are relative to. In DNA this is usually G<->T (6)"),
                            new ElementRule(Parameter.class, true)}
            )
    };

}

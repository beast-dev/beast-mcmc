package dr.evomodelxml.substmodel;

import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.GeneralF81Model;
import dr.xml.*;

import java.util.logging.Logger;

/**
 * Parses a GeneralF81Model
 */
public class GeneralF81ModelParser extends AbstractXMLObjectParser {

    public static final String GENERAL_F81_MODEL = "generalF81Model";
    public static final String FREQUENCIES = "frequencies";

    public String getParserName() {
        return GENERAL_F81_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        FrequencyModel freqModel = null;
        if (xo.hasChildNamed(FREQUENCIES)) {
            XMLObject cxo = xo.getChild(FREQUENCIES);
            freqModel = (FrequencyModel) cxo.getChild(FrequencyModel.class);
        }
        Logger.getLogger("dr.evomodel").info("  General F81 Model from frequencyModel '"+ freqModel.getId() + "' (stateCount=" + freqModel.getFrequencyCount() + ")");

        return new GeneralF81Model(freqModel);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A general F81 model for an arbitrary number of states.";
    }

    public Class getReturnType() {
        return GeneralF81ModelParser.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(FREQUENCIES, FrequencyModel.class)
    };
}

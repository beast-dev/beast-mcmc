package dr.evomodelxml.substmodel;

import dr.evolution.datatype.TwoStateCovarion;
import dr.evomodel.substmodel.AbstractCovarionDNAModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.TwoStateCovarionModel;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * Parses an element from an DOM document into a TwoStateCovarionModel
 */
public class TwoStateCovarionModelParser extends AbstractXMLObjectParser {
    public static final String COVARION_MODEL = "covarionModel";
    public static final String ALPHA = "alpha";
    public static final String SWITCHING_RATE = "switchingRate";

    public String getParserName() {
        return COVARION_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Parameter alphaParameter;
        Parameter switchingRateParameter;

        XMLObject cxo = xo.getChild(AbstractCovarionDNAModel.FREQUENCIES);
        FrequencyModel freqModel = (FrequencyModel) cxo.getChild(FrequencyModel.class);

        TwoStateCovarion dataType = TwoStateCovarion.INSTANCE;  // fancy new datatype courtesy of Helen

        cxo = xo.getChild(ALPHA);
        alphaParameter = (Parameter) cxo.getChild(Parameter.class);

        // alpha must be positive and less than 1.0 because the fast rate is normalized to 1.0
        alphaParameter.addBounds(new Parameter.DefaultBounds(1.0, 0.0, 1));

        cxo = xo.getChild(SWITCHING_RATE);
        switchingRateParameter = (Parameter) cxo.getChild(Parameter.class);

        if (dataType != freqModel.getDataType()) {
            throw new XMLParseException("Data type of " + getParserName() + " element does not match that of its frequencyModel.");
        }

        TwoStateCovarionModel model = new TwoStateCovarionModel(dataType, freqModel, alphaParameter, switchingRateParameter);

        System.out.println(model);

        return model;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A covarion substitution model on binary data and a hidden rate state with two rates.";
    }

    public Class getReturnType() {
        return TwoStateCovarionModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(AbstractCovarionDNAModel.FREQUENCIES, FrequencyModel.class),
            new ElementRule(ALPHA, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class, true)}
            ),
            new ElementRule(SWITCHING_RATE, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class, true)}
            ),
    };
}

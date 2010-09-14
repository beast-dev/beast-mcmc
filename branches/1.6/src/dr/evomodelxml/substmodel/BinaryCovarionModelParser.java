package dr.evomodelxml.substmodel;

import dr.evolution.datatype.TwoStateCovarion;
import dr.evomodel.substmodel.BinaryCovarionModel;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * Parses an element from an DOM document into a TwoStateCovarionModel
 */
public class BinaryCovarionModelParser extends AbstractXMLObjectParser {
    public static final String COVARION_MODEL = "binaryCovarionModel";
    public static final String ALPHA = "alpha";
    public static final String SWITCHING_RATE = "switchingRate";
    public static final String FREQUENCIES = "frequencies";
    public static final String HIDDEN_FREQUENCIES = "hiddenFrequencies"; 

    public String getParserName() {
        return COVARION_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Parameter alphaParameter;
        Parameter switchingRateParameter;

        XMLObject cxo = xo.getChild(FREQUENCIES);
        Parameter frequencies = (Parameter) cxo.getChild(Parameter.class);

        cxo = xo.getChild(HIDDEN_FREQUENCIES);
        Parameter hiddenFrequencies = (Parameter) cxo.getChild(Parameter.class);

        cxo = xo.getChild(ALPHA);
        alphaParameter = (Parameter) cxo.getChild(Parameter.class);

        // alpha must be positive and less than 1.0 because the fast rate is normalized to 1.0
        alphaParameter.addBounds(new Parameter.DefaultBounds(1.0, 0.0, 1));
        hiddenFrequencies.addBounds(new Parameter.DefaultBounds(1.0, 0.0, hiddenFrequencies.getDimension()));
        frequencies.addBounds(new Parameter.DefaultBounds(1.0, 0.0, frequencies.getDimension()));

        cxo = xo.getChild(SWITCHING_RATE);
        switchingRateParameter = (Parameter) cxo.getChild(Parameter.class);

        BinaryCovarionModel model = new BinaryCovarionModel(TwoStateCovarion.INSTANCE,
                frequencies, hiddenFrequencies, alphaParameter, switchingRateParameter);

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
        return BinaryCovarionModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(FREQUENCIES, Parameter.class),
            new ElementRule(HIDDEN_FREQUENCIES, Parameter.class),
            new ElementRule(ALPHA,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class, true)}
            ),
            new ElementRule(SWITCHING_RATE,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class, true)}
            ),
    };


}

package dr.inferencexml.operators;

import dr.xml.*;
import dr.inference.model.Parameter;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.MicrosatelliteAveragingOperator;

/**
 * @author Chieh-Hsi Wu
 *
 * Parser for MicrosatelliteAveragingOperatorParser
 */
public class MicrosatelliteAveragingOperatorParser extends AbstractXMLObjectParser{
    public static final String MODEL_CHOOSE = "modelChoose";
    public static final String DEPENDENCIES = "dependencies";
    public String getParserName() {
        return "msatModelSwitchOperator";
    }
         public Object parseXMLObject(XMLObject xo) throws XMLParseException {
             double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
        Parameter modelChoose = (Parameter) xo.getElementFirstChild(MODEL_CHOOSE);
        Parameter dependencies = (Parameter)xo.getElementFirstChild(DEPENDENCIES);
             return new MicrosatelliteAveragingOperator(modelChoose, dependencies, weight);
    }
         //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************
         public String getParserDescription() {
        return "This element returns a microsatellite averaging operator on a given parameter.";
    }
         public Class getReturnType() {
        return MCMCOperator.class;
    }
         public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }
         private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
            new ElementRule(MODEL_CHOOSE, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(DEPENDENCIES, new XMLSyntaxRule[]{new ElementRule(Parameter.class)})
    };
}

package dr.inferencexml.operators;

import dr.inference.model.Parameter;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.MicrosatelliteModelSelectOperator;
import dr.xml.*;

/**
 * Parser for MicrosatelliteModelSelectOperatorParser
 */
public class MicrosatelliteModelSelectOperatorParser extends AbstractXMLObjectParser {

    public static final String MODEL_INDICATORS = "modelIndicators";
    public static final String MODEL_CHOOSE = "modelChoose";

    public String getParserName() {
        return "msatModelSelectOperator";
    }
         public Object parseXMLObject(XMLObject xo) throws XMLParseException {
             double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
        Parameter modelChoose = (Parameter)xo.getElementFirstChild(MODEL_CHOOSE);
        XMLObject xoInd = xo.getChild(MODEL_INDICATORS);
             int childNum = xoInd.getChildCount();
        System.out.println("There are 12 potential models");
        Parameter[] modelIndicators = new Parameter[childNum];
        for(int i = 0; i < modelIndicators.length; i++){
            modelIndicators[i] = (Parameter)xoInd.getChild(i);
        }
             return new MicrosatelliteModelSelectOperator(modelChoose, modelIndicators, weight);
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
            new ElementRule(MODEL_INDICATORS, new XMLSyntaxRule[]{new ElementRule(Parameter.class,1,Integer.MAX_VALUE)}),
    };
}
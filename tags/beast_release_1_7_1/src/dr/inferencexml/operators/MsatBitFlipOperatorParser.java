package dr.inferencexml.operators;

import dr.xml.*;
import dr.inference.model.Parameter;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.MsatBitFlipOperator;

/**
 * @author Chieh-Hsi Wu
 *
 * Parser for MicrosatelliteAveragingOperatorParser
 */
public class MsatBitFlipOperatorParser extends AbstractXMLObjectParser{
    public static final String MODEL_CHOOSE = "modelChoose";
    public static final String DEPENDENCIES = "dependencies";
    public static final String VARIABLE_INDICES = "variableIndices";

    public String getParserName() {
        return "msatModelSwitchOperator";
    }
         public Object parseXMLObject(XMLObject xo) throws XMLParseException {
             double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
        Parameter modelChoose = (Parameter) xo.getElementFirstChild(MODEL_CHOOSE);
        Parameter dependencies = (Parameter)xo.getElementFirstChild(DEPENDENCIES);
        int[] variableIndices;
            if(xo.hasChildNamed(VARIABLE_INDICES)){

                double[] temp = ((Parameter)xo.getElementFirstChild(VARIABLE_INDICES)).getParameterValues();
                variableIndices = new int[temp.length];
                for(int i = 0; i < temp.length;i++){
                    variableIndices[i] = (int)temp[i];
                }

            }else{
                variableIndices = new int[]{0, 1, 2, 3, 4, 5};
            }

            return new MsatBitFlipOperator(modelChoose, dependencies, weight, variableIndices);
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
            new ElementRule(DEPENDENCIES, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(VARIABLE_INDICES, new XMLSyntaxRule[]{new ElementRule(Parameter.class)},true)
    };
}

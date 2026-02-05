package dr.evomodelxml.epidemiology;
import dr.evomodel.epidemiology.JointCompartmentalModelOperator;
import dr.evomodel.epidemiology.CompartmentalModel;
import dr.inference.model.Parameter;
import dr.inference.operators.SimpleMCMCOperator;
import dr.xml.*;

public class JointCompartmentalModelOperatorParser extends AbstractXMLObjectParser {

    public static final String JOINT_COMPARTMENTAL_MODEL_OPERATOR = "jointCompartmentalModelOperator";
    public static final String OPERATORS = "operators";
    public static final String WEIGHT = "weight";
    public static final String TARGET_ACCEPTANCE = "targetAcceptance";

    public String getParserName() {
        return JOINT_COMPARTMENTAL_MODEL_OPERATOR;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        final double weight = xo.getDoubleAttribute(WEIGHT);

        final double targetProb = xo.getAttribute(TARGET_ACCEPTANCE, 0.2);

        CompartmentalModel compartmentalModel  = (CompartmentalModel) xo.getChild(CompartmentalModel.class);

        if (targetProb <= 0.0 || targetProb >= 1.0)
            throw new RuntimeException("Target acceptance probability must be between 0.0 and 1.0");

        JointCompartmentalModelOperator operator = new JointCompartmentalModelOperator(weight, targetProb, compartmentalModel);

        for(int i = 0; i < xo.getChild(OPERATORS).getChildCount(); i++){
            operator.addOperator((SimpleMCMCOperator) xo.getChild(OPERATORS).getChild(i));
        }

        //for (int i = 0; i < xo.getChildCount(); i++) {
        //    operator.addOperator((SimpleMCMCOperator) xo.getChild(i));
        //}

        return operator;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element represents an arbitrary list of operators for epidemiological parameters, and " +
                "an operator that simulates a compartmental model trajectory given the epidemiological parameters";
    }

    public Class getReturnType() {
        return JointCompartmentalModelOperator.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(CompartmentalModel.class),
            new ElementRule(OPERATORS,
                    new XMLSyntaxRule[]{
                            new ElementRule(SimpleMCMCOperator.class, 1, Integer.MAX_VALUE),}),
            //new ElementRule(SimpleMCMCOperator.class, 1, Integer.MAX_VALUE),
            AttributeRule.newDoubleRule(WEIGHT),
            AttributeRule.newDoubleRule(TARGET_ACCEPTANCE, true)
    };

}

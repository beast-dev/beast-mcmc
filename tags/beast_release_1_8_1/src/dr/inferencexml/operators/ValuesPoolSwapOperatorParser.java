package dr.inferencexml.operators;

import dr.inference.model.ValuesPool;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.ValuesPoolSwapOperator;
import dr.xml.*;

/**
 *
 */
public class ValuesPoolSwapOperatorParser extends AbstractXMLObjectParser {
    public static String VALUESPOOL_OPERATOR = "poolSwapOperator";

    public String getParserName() {
        return VALUESPOOL_OPERATOR;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        final ValuesPool parameter = (ValuesPool) xo.getChild(ValuesPool.class);
        final double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
        final MCMCOperator op = new ValuesPoolSwapOperator(parameter) ;
        op.setWeight(weight);
        return op;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return ".";
    }

    public Class getReturnType() {
        return ValuesPoolSwapOperator.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[] {
                AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
                new ElementRule(ValuesPool.class),
        };
    }

}

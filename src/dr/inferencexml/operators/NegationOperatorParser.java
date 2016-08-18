package dr.inferencexml.operators;

import dr.inference.model.Parameter;
import dr.inference.operators.NegationOperator;
import dr.xml.*;

/**
 * Created by maxryandolinskytolkoff on 8/17/16.
 */
public class NegationOperatorParser extends AbstractXMLObjectParser {
    public final static String NEGATION_OPERATOR = "negationOperator";
    public final static String WEIGHT = "weight";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        Parameter data = (Parameter) xo.getChild(Parameter.class);
        double weight = xo.getDoubleAttribute(WEIGHT);

        return new NegationOperator(data, weight);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[0];
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(Parameter.class),
            AttributeRule.newDoubleRule(WEIGHT),
    };

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return NegationOperatorParser.class;
    }

    @Override
    public String getParserName() {
        return NEGATION_OPERATOR;
    }
}

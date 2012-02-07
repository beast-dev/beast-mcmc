package dr.inferencexml.operators;

import dr.inference.model.Parameter;
import dr.inference.operators.*;
import dr.xml.*;

/**
 */
public class MicrosatUpDownOperatorParser extends AbstractXMLObjectParser {

    public static final String MICROSAT_UP_DOWN_OPERATOR = "microsatUpDownOperator";
    public static final String UP = UpDownOperatorParser.UP;
    public static final String DOWN = UpDownOperatorParser.DOWN;

    public static final String SCALE_FACTOR = ScaleOperatorParser.SCALE_FACTOR;

    public String getParserName() {
        return MICROSAT_UP_DOWN_OPERATOR;
    }

    private Scalable.Default[] getArgs(final XMLObject list) throws XMLParseException {
        Scalable.Default[] args = new Scalable.Default[list.getChildCount()];
        for (int k = 0; k < list.getChildCount(); ++k) {
            final Object child = list.getChild(k);
            if (child instanceof Parameter) {
                args[k] = new Scalable.Default((Parameter) child);
            }

        }
        return args;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        final double scaleFactor = xo.getDoubleAttribute(SCALE_FACTOR);

        final double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);

        final CoercionMode mode = CoercionMode.parseMode(xo);

        final Scalable.Default[] upArgs = getArgs(xo.getChild(UP));
        final Scalable.Default[] dnArgs = getArgs(xo.getChild(DOWN));

        return new MicrosatUpDownOperator(upArgs, dnArgs, scaleFactor, weight, mode);
    }

    public String getParserDescription() {
        return "This element represents an operator that scales two parameters in different directions. " +
                "Each operation involves selecting a scale uniformly at random between scaleFactor and 1/scaleFactor. " +
                "The up parameter is multipled by this scale and the down parameter is divided by this scale.";
    }

    public Class getReturnType() {
        return MicrosatUpDownOperator.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] ee = {
            new ElementRule(Parameter.class, true)
    };

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newDoubleRule(SCALE_FACTOR),
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
            AttributeRule.newBooleanRule(CoercableMCMCOperator.AUTO_OPTIMIZE, true),

            // Allow an arbitrary number of Parameters or Scalables in up or down
            new ElementRule(UP, ee, 1, Integer.MAX_VALUE),
            new ElementRule(DOWN, ee, 1, Integer.MAX_VALUE),
    };
}

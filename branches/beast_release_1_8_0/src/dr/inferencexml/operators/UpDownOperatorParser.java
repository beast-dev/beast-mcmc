package dr.inferencexml.operators;

import dr.inference.model.Parameter;
import dr.inference.operators.*;
import dr.xml.*;

/**
 */
public class UpDownOperatorParser extends AbstractXMLObjectParser {

    public static final String UP_DOWN_OPERATOR = "upDownOperator";
    public static final String UP = "up";
    public static final String DOWN = "down";

    public static final String SCALE_FACTOR = ScaleOperatorParser.SCALE_FACTOR;

    public String getParserName() {
        return UP_DOWN_OPERATOR;
    }

    private Scalable[] getArgs(final XMLObject list) throws XMLParseException {
        Scalable[] args = new Scalable[list.getChildCount()];
        for (int k = 0; k < list.getChildCount(); ++k) {
            final Object child = list.getChild(k);
            if (child instanceof Parameter) {
                args[k] = new Scalable.Default((Parameter) child);
            } else if (child instanceof Scalable) {
                args[k] = (Scalable) child;
            } else {
                XMLObject xo = (XMLObject) child;
                if (xo.hasAttribute("count")) {
                    final int count = xo.getIntegerAttribute("count");

                    final Scalable s = (Scalable) xo.getChild(Scalable.class);
                    args[k] = new Scalable() {

                        public int scale(double factor, int nDims) throws OperatorFailedException {
                            return s.scale(factor, count);
                        }

                        public String getName() {
                            return s.getName() + "(" + count + ")";
                        }
                    };
                } else if (xo.hasAttribute("df")) {
                    final int df = xo.getIntegerAttribute("df");

                    final Scalable s = (Scalable) xo.getChild(Scalable.class);
                    args[k] = new Scalable() {

                        public int scale(double factor, int nDims) throws OperatorFailedException {
                            s.scale(factor, -1);
                            return df;
                        }

                        public String getName() {
                            return s.getName() + "[df=" + df + "]";
                        }
                    };
                }
            }

        }
        return args;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        final double scaleFactor = xo.getDoubleAttribute(SCALE_FACTOR);

        final double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);

        final CoercionMode mode = CoercionMode.parseMode(xo);

        final Scalable[] upArgs = getArgs(xo.getChild(UP));
        final Scalable[] dnArgs = getArgs(xo.getChild(DOWN));

        return new UpDownOperator(upArgs, dnArgs, scaleFactor, weight, mode);
    }

    public String getParserDescription() {
        return "This element represents an operator that scales two parameters in different directions. " +
                "Each operation involves selecting a scale uniformly at random between scaleFactor and 1/scaleFactor. " +
                "The up parameter is multipled by this scale and the down parameter is divided by this scale.";
    }

    public Class getReturnType() {
        return UpDownOperator.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] ee = {
            new ElementRule(Scalable.class, true),
            new ElementRule(Parameter.class, true),
            new ElementRule("scale", new XMLSyntaxRule[]{
                    AttributeRule.newIntegerRule("count", true),
                    AttributeRule.newIntegerRule("df", true),
                    new ElementRule(Scalable.class),
            }, true),
    };

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newDoubleRule(SCALE_FACTOR),
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
            AttributeRule.newBooleanRule(CoercableMCMCOperator.AUTO_OPTIMIZE, true),

            // Allow an arbitrary number of Parameters or Scalable in up or down
            new ElementRule(UP, ee, 1, Integer.MAX_VALUE),
            new ElementRule(DOWN, ee, 1, Integer.MAX_VALUE),
    };

}

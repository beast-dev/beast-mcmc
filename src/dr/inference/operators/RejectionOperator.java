package dr.inference.operators;

import dr.inference.model.Parameter;
import dr.xml.*;

public class RejectionOperator extends SimpleMCMCOperator {

    public interface RejectionProvider extends GibbsOperator {
        double[] getProposedUpdate();

        Parameter getParameter();
    }

    public interface AcceptCondition {
        boolean satisfiesCondition(double[] values);
    }

    public enum SimpleAcceptCondition implements AcceptCondition {
        DescendingAbsoluteValue("descendingAbsoluteValue") {
            @Override
            public boolean satisfiesCondition(double[] values) {
                for (int i = 1; i < values.length; i++) {
                    if (Math.abs(values[i - 1]) < Math.abs(values[i])) {
                        return false;
                    }
                }
                return true;
            }
        },

        DescendingAbsoluteValueSpaced("descendingAbsoluteValueSpaced") {
            @Override
            public boolean satisfiesCondition(double[] values) {
                for (int i = 1; i < values.length; i++) {
                    if (0.9 * Math.abs(values[i - 1]) < Math.abs(values[i])) {
                        return false;
                    }
                }
                return true;
            }
        },

        AlternatingSigns("descendingAlternatingSigns") {
            @Override
            public boolean satisfiesCondition(double[] values) {
                for (int i = 1; i < values.length; i++) {
                    Boolean signa = (values[i] > 0);
                    Boolean signb = (values[i - 1] > 0);
                    if (Math.abs(values[i - 1]) < Math.abs(values[i]) || signa == signb) {
                        return false;
                    }
                }
                return true;
            }
        };

        private final String name;

        SimpleAcceptCondition(String name) {
            this.name = name;
        }

        public String getName() {
            return name;
        }

        @Override
        public abstract boolean satisfiesCondition(double[] values);
    }

    private final RejectionProvider gibbsOp;
    private final AcceptCondition condition;

    public RejectionOperator(RejectionProvider gibbsOp, AcceptCondition condition, double weight) {
        setWeight(weight);
        this.gibbsOp = gibbsOp;
        this.condition = condition;
    }

    @Override
    public String getOperatorName() {
        return REJECTION_OPERATOR;
    }

    @Override
    public double doOperation() {
        double[] values = gibbsOp.getProposedUpdate();
        if (condition.satisfiesCondition(values)) {

            Parameter parameter = gibbsOp.getParameter();
            for (int i = 0; i < parameter.getDimension(); i++) {
                parameter.setParameterValueQuietly(i, values[i]);
            }

            parameter.fireParameterChangedEvent();

            return Double.POSITIVE_INFINITY;
        }
        return Double.NEGATIVE_INFINITY;
    }


    private static final String REJECTION_OPERATOR = "rejectionOperator";
    private static final String CONDITION = "condition";

    public static AbstractXMLObjectParser PARSER = new AbstractXMLObjectParser() {
        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            RejectionProvider gibbsOp = (RejectionProvider) xo.getChild(RejectionProvider.class);
            double weight = xo.getAttribute(WEIGHT, gibbsOp.getWeight());

            AcceptCondition condition = null;
            String stringCondition = xo.getStringAttribute(CONDITION);
            if (stringCondition == null) {
                condition = (AcceptCondition) xo.getChild(AcceptCondition.class);
            } else {
                for (SimpleAcceptCondition simpleCondition : SimpleAcceptCondition.values()) {
                    if (stringCondition.equalsIgnoreCase(simpleCondition.getName())) {
                        condition = simpleCondition;
                        break;
                    }
                }
                if (condition == null) {
                    throw new XMLParseException("Unrecognized condition type: " + stringCondition);
                }
            }
            return new RejectionOperator(gibbsOp, condition, weight);
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return new XMLSyntaxRule[]{
                    new ElementRule(RejectionProvider.class),
                    new XORRule(
                            AttributeRule.newStringRule(CONDITION),
                            new ElementRule(AcceptCondition.class)
                    ),
                    AttributeRule.newDoubleRule(WEIGHT, true)
            };
        }

        @Override
        public String getParserDescription() {
            return "A rejection sampler that always accepts a Gibbs proposal if a condition is met (and always rejects otherwise)";
        }

        @Override
        public Class getReturnType() {
            return RejectionOperator.class;
        }

        @Override
        public String getParserName() {
            return REJECTION_OPERATOR;
        }
    };


}

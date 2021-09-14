package dr.inference.operators.rejection;

import dr.xml.*;

public class DescendingAndSpacedCondition implements AcceptCondition {
    private final double spacing;

    DescendingAndSpacedCondition(double spacing) {
        this.spacing = spacing;
    }

    @Override
    public boolean satisfiesCondition(double[] values) {
        for (int i = 1; i < values.length; i++) {
            if (spacing * Math.abs(values[i - 1]) < Math.abs(values[i])) {
                return false;
            }
        }
        return true;
    }

    private static final String DESCENDING_AND_SPACED = "descendingAndSpaced";
    private static final String SPACING = "spacing";

    public static AbstractXMLObjectParser PARSER = new AbstractXMLObjectParser() {
        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            double spacing = xo.getDoubleAttribute(SPACING);
            if (spacing < 0.0 || spacing > 1.0) {
                throw new XMLParseException("Attribute '" + SPACING + "' must be between 0 and 1.");
            }

            return new DescendingAndSpacedCondition(spacing);
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return new XMLSyntaxRule[]{
                    AttributeRule.newDoubleRule(SPACING)
            };
        }

        @Override
        public String getParserDescription() {
            return "Condition requiring parameter to have descending absolute values with some minimum spacing.";
        }

        @Override
        public Class getReturnType() {
            return DescendingAndSpacedCondition.class;
        }

        @Override
        public String getParserName() {
            return DESCENDING_AND_SPACED;
        }
    };
}


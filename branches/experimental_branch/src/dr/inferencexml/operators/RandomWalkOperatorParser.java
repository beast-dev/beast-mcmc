package dr.inferencexml.operators;

import dr.inference.model.Parameter;
import dr.inference.operators.CoercableMCMCOperator;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.RandomWalkOperator;
import dr.xml.*;

/**
 */
public class RandomWalkOperatorParser extends AbstractXMLObjectParser {

    public static final String RANDOM_WALK_OPERATOR = "randomWalkOperator";
    public static final String WINDOW_SIZE = "windowSize";
    public static final String UPDATE_INDEX = "updateIndex";
    public static final String UPPER = "upper";
    public static final String LOWER = "lower";

    public static final String BOUNDARY_CONDITION = "boundaryCondition";

        public String getParserName() {
            return RANDOM_WALK_OPERATOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            CoercionMode mode = CoercionMode.parseMode(xo);

            double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
            double windowSize = xo.getDoubleAttribute(WINDOW_SIZE);
            Parameter parameter = (Parameter) xo.getChild(Parameter.class);

            Double lower = null;
            Double upper = null;

            if (xo.hasAttribute(LOWER)) {
                lower = xo.getDoubleAttribute(LOWER);
            }

            if (xo.hasAttribute(UPPER)) {
                upper = xo.getDoubleAttribute(UPPER);
            }

            RandomWalkOperator.BoundaryCondition condition = RandomWalkOperator.BoundaryCondition.valueOf(
                    xo.getAttribute(BOUNDARY_CONDITION, RandomWalkOperator.BoundaryCondition.reflecting.name()));

            if (xo.hasChildNamed(UPDATE_INDEX)) {
                XMLObject cxo = xo.getChild(UPDATE_INDEX);
                Parameter updateIndex = (Parameter) cxo.getChild(Parameter.class);
                if (updateIndex.getDimension() != parameter.getDimension())
                    throw new RuntimeException("Parameter to update and missing indices must have the same dimension");
                return new RandomWalkOperator(parameter, updateIndex, windowSize, condition,
                        weight, mode, lower, upper);
            }

            return new RandomWalkOperator(parameter, null, windowSize, condition, weight, mode, lower, upper);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element returns a random walk operator on a given parameter.";
        }

        public Class getReturnType() {
            return MCMCOperator.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newDoubleRule(WINDOW_SIZE),
                AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
                AttributeRule.newDoubleRule(LOWER, true),
                AttributeRule.newDoubleRule(UPPER, true),
                AttributeRule.newBooleanRule(CoercableMCMCOperator.AUTO_OPTIMIZE, true),
                new ElementRule(UPDATE_INDEX,
                        new XMLSyntaxRule[] {
                                new ElementRule(Parameter.class),
                        },true),
                new StringAttributeRule(BOUNDARY_CONDITION, null, RandomWalkOperator.BoundaryCondition.values(), true),
                new ElementRule(Parameter.class)
        };
}

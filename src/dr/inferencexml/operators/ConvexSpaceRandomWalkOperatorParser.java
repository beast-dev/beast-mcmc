package dr.inferencexml.operators;

import dr.inference.model.Parameter;
import dr.inference.operators.ConvexSpaceRandomWalkOperator;
import dr.inference.operators.MCMCOperator;
import dr.math.distributions.ConvexSpaceRandomGenerator;
import dr.xml.*;

public class ConvexSpaceRandomWalkOperatorParser extends AbstractXMLObjectParser {


    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        Parameter parameter = (Parameter) xo.getChild(Parameter.class);
        ConvexSpaceRandomGenerator generator =
                (ConvexSpaceRandomGenerator) xo.getChild(ConvexSpaceRandomGenerator.class);

        if (!generator.isUniform()) {
            throw new XMLParseException("sample distribution must be uniform over its support");
        }

        double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);

        double windowSize = xo.getAttribute(ConvexSpaceRandomWalkOperator.WINDOW_SIZE, 1.0);
        if (windowSize > 1.0) {
            throw new XMLParseException(ConvexSpaceRandomWalkOperator.WINDOW_SIZE + " must be between 0 and 1");
        }

        return new ConvexSpaceRandomWalkOperator(parameter, generator, windowSize, weight);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                new ElementRule(Parameter.class),
                new ElementRule(ConvexSpaceRandomGenerator.class),
                AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
                AttributeRule.newDoubleRule(ConvexSpaceRandomWalkOperator.WINDOW_SIZE, true)
        };
    }

    @Override
    public String getParserDescription() {
        return "operator that first samples uniformly from some space then updates the parameter to a point along" +
                " the line from its current value to the sampled one";
    }

    @Override
    public Class getReturnType() {
        return ConvexSpaceRandomWalkOperator.class;
    }

    @Override
    public String getParserName() {
        return ConvexSpaceRandomWalkOperator.CONVEX_RW;
    }
}

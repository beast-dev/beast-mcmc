package dr.inferencexml.operators;

import dr.inference.model.BoundedSpace;
import dr.inference.model.Parameter;
import dr.inference.operators.ConvexSpaceRandomWalkOperator;
import dr.inference.operators.MCMCOperator;
import dr.xml.*;

public class ConvexSpaceRandomWalkOperatorParser extends AbstractXMLObjectParser {


    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        Parameter parameter = (Parameter) xo.getChild(Parameter.class);
        BoundedSpace space =
                (BoundedSpace) xo.getChild(BoundedSpace.class);


        double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);

        double windowSize = xo.getAttribute(ConvexSpaceRandomWalkOperator.WINDOW_SIZE, 1.0);
        if (windowSize > 1.0) {
            throw new XMLParseException(ConvexSpaceRandomWalkOperator.WINDOW_SIZE + " must be between 0 and 1");
        }

        final Parameter updateIndex;

        if (xo.hasChildNamed(RandomWalkOperatorParser.UPDATE_INDEX)) {
            XMLObject cxo = xo.getChild(RandomWalkOperatorParser.UPDATE_INDEX);
            updateIndex = (Parameter) cxo.getChild(Parameter.class);
        } else {
            updateIndex = null;
        }

        return new ConvexSpaceRandomWalkOperator(parameter, space, updateIndex, windowSize, weight);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                new ElementRule(Parameter.class),
                new ElementRule(BoundedSpace.class),
                new ElementRule(RandomWalkOperatorParser.UPDATE_INDEX,
                        new XMLSyntaxRule[]{
                                new ElementRule(Parameter.class)
                        }, true
                ),
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

package dr.inferencexml.distribution;

import dr.inference.distribution.RandomWalkGenerator;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * Created by mkarcher on 4/5/17.
 */
public class RandomWalkGeneratorParser extends AbstractXMLObjectParser {
    public static final String RANDOM_WALK_GENERATOR = "randomWalkGenerator";
    public static final String FIRST_ELEM_PREC = "firstElementPrecision";
    public static final String PREC = "precision";
    public static final String DIM = "dimension";

    @Override
    public String getParserName() { return RANDOM_WALK_GENERATOR; }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        XMLObject cxo = xo.getChild(FIRST_ELEM_PREC);
        Parameter firstElementPrecision = (Parameter) cxo.getChild(Parameter.class);

        cxo = xo.getChild(PREC);
        Parameter prec = (Parameter) cxo.getChild(Parameter.class);

        cxo = xo.getChild(DIM); // May need to adapt to multiple trees, a la CoalescentLikelihoodParser
        Integer dim = (Integer) cxo.getChild(Integer.class);

        return new RandomWalkGenerator(dim.intValue(), firstElementPrecision, prec);
    }

    @Override
    public String getParserDescription() {
        return "This element generates a regular Gaussian random walk.";
    }

    @Override
    public Class getReturnType() {
        return RandomWalkGenerator.class;
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(FIRST_ELEM_PREC, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }, "The precision for the first element of the regular random walk"),

            new ElementRule(PREC, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }, "The precision for the relationship between adjacent elements in the random walk"),

            new ElementRule(DIM, new XMLSyntaxRule[]{
                    new ElementRule(Integer.class)
            }, "The number of elements in the random walk")
    };
}

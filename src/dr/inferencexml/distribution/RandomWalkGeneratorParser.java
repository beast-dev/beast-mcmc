package dr.inferencexml.distribution;

import dr.inference.distribution.RandomWalkGenerator;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * Created by mkarcher on 4/5/17.
 */
public class RandomWalkGeneratorParser extends AbstractXMLObjectParser {
    public static final String RANDOM_WALK_GENERATOR = "randomWalkGenerator";
    public static final String DATA = "data";
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

        cxo = xo.getChild(DATA); // May need to adapt to multiple trees, a la CoalescentLikelihoodParser
        Parameter data = (Parameter) cxo.getChild(Parameter.class);

        return new RandomWalkGenerator(data, firstElementPrecision, prec);
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
            new ElementRule(DATA, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }, "The data to evaluate the density"),


            new ElementRule(FIRST_ELEM_PREC, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }, "The precision for the first element of the regular random walk"),

            new ElementRule(PREC, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }, "The precision for the relationship between adjacent elements in the random walk"),

//            AttributeRule.newIntegerRule(DIM)
    };
}

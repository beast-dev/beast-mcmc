package dr.inferencexml.operators;

import dr.inference.model.Parameter;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.OneZeroOneShuffleOperator;
import dr.xml.*;

/**
 * <oneZeroOneShuffleOperator id="..." weight="..." tol="...">
 *   <alphaRates><parameter idref="ind.rewardRates"/></alphaRates>
 *   <extremeIndex><parameter idref="extremeIndex"/></extremeIndex>
 * </oneZeroOneShuffleOperator>
 */

/*
 * @author: Filippo Monti
 */
public final class OneZeroOneShuffleOperatorParser extends AbstractXMLObjectParser {

    public static final String OPERATOR = "oneZeroOneShuffleOperator";

    private static final String WEIGHT = "weight";
    private static final String TOL = "tol";

    private static final String ALPHA_RATES = "alphaRates";
    private static final String EXTREME_INDEX = "extremeIndex";

    @Override
    public String getParserName() {
        return OPERATOR;
    }

    @Override
    public Object parseXMLObject(final XMLObject xo) throws XMLParseException {

        final double weight = xo.getDoubleAttribute(WEIGHT);
        final double tol = xo.hasAttribute(TOL) ? xo.getDoubleAttribute(TOL) : 0.0;

        final Parameter alphaRates = (Parameter) xo.getElementFirstChild(ALPHA_RATES);
        final Parameter extremeIndex = (Parameter) xo.getElementFirstChild(EXTREME_INDEX);


            final OneZeroOneShuffleOperator op =
                    new OneZeroOneShuffleOperator(alphaRates, extremeIndex, weight, tol);

            return op;

    }

    @Override
    public String getParserDescription() {
        return "Relabels which alphaRates entries are forced to be 0 and 1 via a length-2 extremeIndex parameter.";
    }

    @Override
    public Class getReturnType() {
        return MCMCOperator.class;
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private static final XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
            AttributeRule.newDoubleRule(WEIGHT),
            AttributeRule.newDoubleRule(TOL, true),

            new ElementRule(ALPHA_RATES,
                    new XMLSyntaxRule[] { new ElementRule(Parameter.class) },
                    "The K-length alphaRates parameter (e.g., ind.rewardRates).", false),

            new ElementRule(EXTREME_INDEX,
                    new XMLSyntaxRule[] { new ElementRule(Parameter.class) },
                    "Length-2 integer-like parameter giving indices of the 0 and 1 labeled extremes.", false),
    };
}
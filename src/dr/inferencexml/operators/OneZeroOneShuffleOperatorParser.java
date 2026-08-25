package dr.inferencexml.operators;

import dr.evomodel.branchratemodel.RewardRates;
import dr.inference.operators.MCMCOperator;
import dr.xml.*;

/**
 * Disabled: <oneZeroOneShuffleOperator> has a reachability bug (0/1998
 * acceptance observed in practice, see ctmc_bm4d_timeseries scenario a) and
 * must not be used. Retained only so any XML still referencing it fails
 * fast with a clear message instead of silently proposing an unreachable
 * move. Use <oneZeroOneShuffleGibbsOperator> instead.
 */

/*
 * @author: Filippo Monti
 */
public final class OneZeroOneShuffleOperatorParser extends AbstractXMLObjectParser {

    public static final String OPERATOR = "oneZeroOneShuffleOperator";

    private static final String WEIGHT = "weight";
    private static final String TOL = "tol";

    @Override
    public String getParserName() {
        return OPERATOR;
    }

    @Override
    public Object parseXMLObject(final XMLObject xo) throws XMLParseException {
        throw new XMLParseException(
                "<" + OPERATOR + "> has a reachability bug (0/1998 acceptance observed in " +
                        "practice, see ctmc_bm4d_timeseries scenario a) and must not be used. " +
                        "Use <" + OneZeroOneShuffleGibbsOperatorParser.OPERATOR + "> instead."
        );
    }


    @Override
    public String getParserDescription() {
        return "Disabled (reachability bug) -- use oneZeroOneShuffleGibbsOperator instead.";
    }

    @Override
    public Class getReturnType() {
        return MCMCOperator.class;
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private static final XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            AttributeRule.newDoubleRule(WEIGHT),
            AttributeRule.newDoubleRule(TOL, true),

            new ElementRule(RewardRates.class)
    };
}

package dr.inferencexml.operators;

import dr.evomodel.branchratemodel.RewardRates;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.OneZeroOneShuffleOperator;
import dr.xml.*;

/**
 * <oneZeroOneShuffleOperator id="..." weight="..." tol="...">
 *   <rewardRates>
 *     <fixedValues>
 *       <parameter idref="ind.rewardRates.fixedValues"/>
 *     </fixedValues>
 *     <varyingValues>
 *       <parameter idref="ind.rewardRates.varyingValues"/>
 *     </varyingValues>
 *     <stateIndices>
 *       <parameter idref="ind.rewardRates.stateIndices"/>
 *     </stateIndices>
 *   </rewardRates>
 * </oneZeroOneShuffleOperator>
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

        final double weight = xo.getDoubleAttribute(WEIGHT);
        final double tol = xo.hasAttribute(TOL) ? xo.getDoubleAttribute(TOL) : 0.0;

        final RewardRates rewardRates = (RewardRates) xo.getChild(RewardRates.class);

        final OneZeroOneShuffleOperator op =
                new OneZeroOneShuffleOperator(
                        rewardRates.getValues(),
                        rewardRates.getStateIndices(),
                        weight,
                        tol);

        return op;

    }


    @Override
    public String getParserDescription() {
        return "Swaps one fixed extreme reward-rate label (0 or 1) with one non-extreme reward-rate label in the state-to-reward-rate mapping.";
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

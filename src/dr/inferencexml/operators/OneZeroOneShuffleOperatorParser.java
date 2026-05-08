package dr.inferencexml.operators;

import dr.evomodelxml.branchratemodel.RewardsAwareMixtureBranchRatesParser;
import dr.inference.model.Parameter;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.OneZeroOneShuffleOperator;
import dr.xml.*;

/**
 * <oneZeroOneShuffleOperator id="..." weight="..." tol="...">
 *   <rewardRates>
 *     <values>
 *       <parameter idref="ind.rewardRates.values"/>
 *     </values>
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

        final XMLObject rewardRatesObj = xo.getChild(RewardsAwareMixtureBranchRatesParser.REWARD_RATES);
        final Parameter rewardRatesValues = (Parameter)
                rewardRatesObj.getElementFirstChild(RewardsAwareMixtureBranchRatesParser.REWARD_RATES_VALUES);
        final Parameter rewardRatesMapping = (Parameter)
                rewardRatesObj.getElementFirstChild(RewardsAwareMixtureBranchRatesParser.REWARD_RATES_MAPPING);

        final OneZeroOneShuffleOperator op =
                new OneZeroOneShuffleOperator(rewardRatesValues, rewardRatesMapping, weight, tol);

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

            new ElementRule(RewardsAwareMixtureBranchRatesParser.REWARD_RATES, new XMLSyntaxRule[]{
                    new ElementRule(RewardsAwareMixtureBranchRatesParser.REWARD_RATES_VALUES, new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class)
                    }),
                    new ElementRule(RewardsAwareMixtureBranchRatesParser.REWARD_RATES_MAPPING, new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class)
                    })
            })
    };
}

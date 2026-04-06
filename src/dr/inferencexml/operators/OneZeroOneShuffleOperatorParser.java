package dr.inferencexml.operators;

import dr.evomodelxml.branchratemodel.RewardsAwareMixtureBranchRatesParser;
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

    private static final XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            AttributeRule.newDoubleRule(WEIGHT),
            AttributeRule.newDoubleRule(TOL, true),

//            new ElementRule(ALPHA_RATES,
//                    new XMLSyntaxRule[] { new ElementRule(Parameter.class) },
//                    "The K-length alphaRates parameter (e.g., ind.rewardRates).", false),

            new ElementRule(RewardsAwareMixtureBranchRatesParser.REWARD_RATES, new XMLSyntaxRule[]{
                    new ElementRule(RewardsAwareMixtureBranchRatesParser.REWARD_RATES_VALUES, new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class)
                    }),
//                        new ElementRule(RewardsAwareMixtureBranchRatesParser.REWARD_RATES_VARYING_VALUES, new XMLSyntaxRule[] {
//                                new ElementRule(Parameter.class)
//                        }),
                    new ElementRule(RewardsAwareMixtureBranchRatesParser.REWARD_RATES_MAPPING, new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class)
                    })
            })
    };
}
package dr.evomodelxml.branchratemodel;

import dr.evomodel.branchratemodel.RewardRates;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Parameter;
import dr.xml.*;

/*
 * @author Filippo Monti
 */
public final class RewardRatesParser extends AbstractXMLObjectParser {

    public static final String REWARD_RATES = "rewardRates";
    public static final String FIXED_VALUES = "fixedValues";
    public static final String VARYING_VALUES = "varyingValues";
    public static final String STATE_INDICES = "stateIndices";

    @Override
    public String getParserName() {
        return REWARD_RATES;
    }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        final Parameter fixedValues = parseParameter(xo, FIXED_VALUES);
        final Parameter varyingValues = parseParameter(xo, VARYING_VALUES);
        final Parameter values = buildValuesParameter(xo, fixedValues, varyingValues);
        final Parameter stateIndices = parseParameter(xo, STATE_INDICES);

        validateValues(values);
        validatePermutation(stateIndices);
        if (values.getDimension() != stateIndices.getDimension()) {
            throw new XMLParseException("rewardRates values dim must match rewardRates stateIndices dim");
        }
        if (varyingValues.getDimension() != values.getDimension() - 2) {
            throw new XMLParseException("rewardRates varyingValues dim must match values dim minus 2");
        }

        return new RewardRates(values, fixedValues, varyingValues, stateIndices);
    }

    private static Parameter parseParameter(XMLObject xo, String name) throws XMLParseException {
        final Object child = xo.getElementFirstChild(name);
        if (!(child instanceof Parameter)) {
            throw new XMLParseException("rewardRates/" + name + " must contain a parameter.");
        }
        return (Parameter) child;
    }

    private static Parameter buildValuesParameter(XMLObject xo, Parameter fixedValues, Parameter varyingValues)
            throws XMLParseException {
        if (fixedValues.getDimension() != 2) {
            throw new XMLParseException("rewardRates/fixedValues must have dimension 2.");
        }
        final CompoundParameter values = new CompoundParameter(
                xo.hasId() ? xo.getId() + ".values" : null);
        values.addParameter(fixedValues);
        values.addParameter(varyingValues);
        return values;
    }

    public static void validatePermutation(Parameter mapping) throws XMLParseException {
        final int K = mapping.getDimension();
        final boolean[] seen = new boolean[K];

        for (int i = 0; i < K; i++) {
            final double v = mapping.getParameterValue(i);
            final int idx = (int) Math.round(v);
            if (Math.abs(v - idx) > 1e-10) {
                throw new XMLParseException("rewardRates stateIndices must contain integer values.");
            }
            if (idx < 0 || idx >= K) {
                throw new XMLParseException("rewardRates stateIndices entries must be in [0, " + (K - 1) + "].");
            }
            if (seen[idx]) {
                throw new XMLParseException("rewardRates stateIndices contains duplicate value: " + idx);
            }
            seen[idx] = true;
        }

        for (int i = 0; i < K; i++) {
            if (!seen[i]) {
                throw new XMLParseException("rewardRates stateIndices is missing value: " + i);
            }
        }
    }

    public static void validateValues(Parameter values) throws XMLParseException {
        final int K = values.getDimension();
        if (K < 2) {
            throw new XMLParseException("rewardRates values must have dimension at least 2 (for 0 and 1).");
        }

        final double tol = 1e-10;
        final double v0 = values.getParameterValue(0);
        final double v1 = values.getParameterValue(1);

        if (Math.abs(v0 - 0.0) > tol) {
            throw new XMLParseException("rewardRates values[0] must be 0, found: " + v0);
        }
        if (Math.abs(v1 - 1.0) > tol) {
            throw new XMLParseException("rewardRates values[1] must be 1, found: " + v1);
        }
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private static final XMLSyntaxRule[] rules = {
            new ElementRule(FIXED_VALUES, new XMLSyntaxRule[] {
                    new ElementRule(Parameter.class)
            }),
            new ElementRule(VARYING_VALUES, new XMLSyntaxRule[] {
                    new ElementRule(Parameter.class)
            }),
            new ElementRule(STATE_INDICES, new XMLSyntaxRule[] {
                    new ElementRule(Parameter.class)
            })
    };

    @Override
    public String getParserDescription() {
        return "Reward-rate values, free varying values, and state-to-reward-rate mapping.";
    }

    @Override
    public Class getReturnType() {
        return RewardRates.class;
    }
}

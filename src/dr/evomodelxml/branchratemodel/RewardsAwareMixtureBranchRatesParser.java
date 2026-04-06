package dr.evomodelxml.branchratemodel;

import dr.evomodel.branchratemodel.ArbitraryBranchRates;
import dr.evomodel.branchratemodel.RewardsAwareMixtureBranchRates;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.tree.TreeParameterModel;
import dr.inference.model.Parameter;
import dr.xml.*;

/*
 * @author Filippo Monti
 */
public final class RewardsAwareMixtureBranchRatesParser extends AbstractXMLObjectParser {

    public static final String NAME = "rewardsAwareMixtureBranchRates";
    public static final String ATOMIC = "atomic";
    public static final String CTS = "ctsParameter";
    public static final String INDICATOR = "indicator";
    public static final String INCLUDE_ROOT = "includeRoot";
    public static final String ATOMS_INDICES = "atomIndices";
    public static final String REWARD_RATES = "rewardRates";
    public static final String REWARD_RATES_VALUES = "values";
    public static final String REWARD_RATES_VARYING_VALUES = "varyingValues";
    public static final String REWARD_RATES_MAPPING = "stateIndices";

    @Override
    public String getParserName() { return NAME; }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);


        Parameter indicator = (Parameter) xo.getChild(INDICATOR).getChild(Parameter.class);
        Parameter cts = (Parameter) xo.getChild(CTS).getChild(Parameter.class);

        XMLObject atomicObj = xo.getChild(ATOMIC);

        Parameter atomIndices = (Parameter) atomicObj.getElementFirstChild(ATOMS_INDICES);

        XMLObject rewardsRatesObj = atomicObj.getChild(REWARD_RATES);
        Parameter rewardRatesValues = (Parameter) rewardsRatesObj.getElementFirstChild(REWARD_RATES_VALUES);
        validateRewardRatesValues(rewardRatesValues);
        Parameter rewardRatesVaryingValues = (Parameter) rewardsRatesObj.getElementFirstChild(REWARD_RATES_VARYING_VALUES);
        Parameter rewardRatesMapping = (Parameter) rewardsRatesObj.getElementFirstChild(REWARD_RATES_MAPPING);
        validatePermutation(rewardRatesMapping);
        if (rewardRatesValues.getDimension() != rewardRatesMapping.getDimension()) {
            throw new XMLParseException("rewardRates values dim must match the number of states indicated by rewardRates mapping dim");
        }
        if (rewardRatesVaryingValues.getDimension() != rewardRatesValues.getDimension() - 2) {
            throw new XMLParseException("rewardRates varying values dim must match the number of states minus 2");
        }

        TreeParameterModel.Type includeRoot = xo.getAttribute(INCLUDE_ROOT, false)
                ? TreeParameterModel.Type.WITH_ROOT
                : TreeParameterModel.Type.WITHOUT_ROOT;

        final int nBranches = tree.getNodeCount() - 1;
        if (includeRoot == TreeParameterModel.Type.WITHOUT_ROOT) {
            if (cts.getDimension() != nBranches) throw new XMLParseException("cts dim must be nodeCount-1");
            if (indicator.getDimension() != nBranches) throw new XMLParseException("indicator dim must be nodeCount-1");
            if (atomIndices.getDimension() != nBranches) throw new XMLParseException("atomic dim must be nodeCount-1");
        } else {
            if (cts.getDimension() != nBranches + 1) throw new XMLParseException("cts dim must be nodeCount");
            if (indicator.getDimension() != nBranches + 1) throw new XMLParseException("indicator dim must be nodeCount");
            if (atomIndices.getDimension() != nBranches + 1) throw new XMLParseException("atomic dim must be nodeCount");
        }

        final ArbitraryBranchRates.BranchRateTransform transform = ArbitraryBranchRatesParser.parseTransform(xo);

        return new RewardsAwareMixtureBranchRates(tree, cts, indicator, atomIndices,
                rewardRatesValues, rewardRatesVaryingValues, rewardRatesMapping,
                transform, false, includeRoot);
    }

    public static void validatePermutation(Parameter mapping) throws XMLParseException {
        final int K = mapping.getDimension();
        final boolean[] seen = new boolean[K];

        for (int i = 0; i < K; i++) {
            double v = mapping.getParameterValue(i);

            // check integer-like
            int idx = (int) Math.round(v);
            if (Math.abs(v - idx) > 1e-10) {
                throw new XMLParseException("rewardRatesMapping must contain integer values.");
            }

            // check bounds
            if (idx < 0 || idx >= K) {
                throw new XMLParseException("rewardRatesMapping entries must be in [0, " + (K - 1) + "].");
            }

            // check uniqueness
            if (seen[idx]) {
                throw new XMLParseException("rewardRatesMapping contains duplicate value: " + idx);
            }

            seen[idx] = true;
        }

        // optional: this is redundant but makes intent explicit
        for (int i = 0; i < K; i++) {
            if (!seen[i]) {
                throw new XMLParseException("rewardRatesMapping is missing value: " + i);
            }
        }
    }

    private static void validateRewardRatesValues(Parameter values) throws XMLParseException {
        final int K = values.getDimension();

        if (K < 2) {
            throw new XMLParseException("rewardRatesValues must have dimension at least 2 (for 0 and 1).");
        }

        final double tol = 1e-10;

        double v0 = values.getParameterValue(0);
        double v1 = values.getParameterValue(1);

        if (Math.abs(v0 - 0.0) > tol) {
            throw new XMLParseException("rewardRatesValues[0] must be 0, found: " + v0);
        }

        if (Math.abs(v1 - 1.0) > tol) {
            throw new XMLParseException("rewardRatesValues[1] must be 1, found: " + v1);
        }
    }

    @Override
    public String getParserDescription() {
        return "Mixture branch-rate model: per-branch continuous rate or atomic rate selected via an IndexedParameter, controlled by an indicator.";
    }

    @Override
    public Class getReturnType() { return RewardsAwareMixtureBranchRates.class; }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[] {
                new ElementRule(TreeModel.class),

                new ElementRule(CTS, new XMLSyntaxRule[] {
                        new ElementRule(Parameter.class)
                }),

                new ElementRule(INDICATOR, new XMLSyntaxRule[] {
                        new ElementRule(Parameter.class)
                }),

                new ElementRule(ATOMIC, new XMLSyntaxRule[] {
                        new ElementRule(ATOMS_INDICES, new XMLSyntaxRule[] {
                                new ElementRule(Parameter.class)
                        }),
                        new ElementRule(REWARD_RATES, new XMLSyntaxRule[] {
                                new ElementRule(REWARD_RATES_VALUES, new XMLSyntaxRule[] {
                                        new ElementRule(Parameter.class)
                                }),
                                new ElementRule(REWARD_RATES_VARYING_VALUES, new XMLSyntaxRule[] {
                                        new ElementRule(Parameter.class)
                                }),
                                new ElementRule(REWARD_RATES_MAPPING, new XMLSyntaxRule[] {
                                        new ElementRule(Parameter.class)
                                })
                        })
                }),

                AttributeRule.newBooleanRule(INCLUDE_ROOT, true),

                new ElementRule(ArbitraryBranchRatesParser.SCALE,
                        Parameter.class, "optional scale parameter", true),

                new ElementRule(ArbitraryBranchRatesParser.LOCATION,
                        Parameter.class, "optional location parameter", true),

                new ElementRule(ArbitraryBranchRatesParser.RANDOM_INDICATOR,
                        new XMLSyntaxRule[] {
                                new ElementRule(Parameter.class),
                        }, true),

                AttributeRule.newBooleanRule(ArbitraryBranchRatesParser.SHRINKAGE, true),
        };
    }

}
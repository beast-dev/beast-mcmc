package dr.evomodel.branchratemodel;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.tree.TreeParameterModel;
import dr.inference.model.IndexedParameter;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

/**
 * Mixture branch-rate model:
 *
 *   rho_b = atomic_b  if indicator_b == 1
 *         = cts_b     if indicator_b == 0
 *
 * Continuous rates and transform behavior are inherited from ArbitraryBranchRates.
 * Only the raw branch-rate selection and the derivative masking are customized.
 */
/*
 * @author Filippo Monti
 */
public final class RewardsAwareMixtureBranchRates extends ArbitraryBranchRates {

    public static final String ID = "rewardsAwareMixtureBranchRates";

    private final Parameter indicator;
    private final Parameter atomIndices;
    private final Parameter rewardRatesValues;
    private final Parameter rewardRatesMapping;
    private final Parameter rewardRatesInternal;


    public RewardsAwareMixtureBranchRates(TreeModel tree,
                                          Parameter ctsParameter,
                                          Parameter indicator,
                                          Parameter atomIndices,
                                          Parameter rewardRatesValues,
                                          Parameter rewardRatesInternal,
                                          Parameter rewardRatesMapping,
                                          BranchRateTransform transform,
                                          boolean setRates,
                                          TreeParameterModel.Type includeRoot) {
        super(ID, tree, ctsParameter,
                transform == null ? new BranchRateTransform.None() : transform,
                setRates, includeRoot);

        if (indicator == null) {
            throw new IllegalArgumentException("indicator must be non-null");
        }
        if (atomIndices == null) {
            throw new IllegalArgumentException("atomIndices must be non-null");
        }
        if (rewardRatesValues == null) {
            throw new IllegalArgumentException("rewardRatesValues must be non-null");
        }
        if (rewardRatesMapping == null) {
            throw new IllegalArgumentException("rewardRatesMapping must be non-null");
        }
        if (rewardRatesInternal == null) {
            throw new IllegalArgumentException("rewardRatesInternal must be non-null");
        }

        this.indicator = indicator;
        this.atomIndices = atomIndices;
        this.rewardRatesValues = rewardRatesValues;
        this.rewardRatesInternal = rewardRatesInternal;
        this.rewardRatesMapping = rewardRatesMapping;

        final int expected = ctsParameter.getDimension();

        if (indicator.getDimension() != expected) {
            throw new IllegalArgumentException(
                    "indicator dim must match ctsParameter dim (" + expected + ") but is " +
                            indicator.getDimension());
        }
        if (atomIndices.getDimension() != expected) {
            throw new IllegalArgumentException(
                    "atomIndices dim must match ctsParameter dim (" + expected + ") but is " +
                            atomIndices.getDimension());
        }

        addVariable(indicator);
        addVariable(atomIndices);
        addVariable(rewardRatesValues);
        addVariable(rewardRatesInternal);
        addVariable(rewardRatesMapping);
    }

    @Override
    public double getUntransformedBranchRate(final Tree tree, final NodeRef node) {
        final int p = getParameterIndexFromNode(node);

        if (isOne(indicator.getParameterValue(p))) {
            final int stateIndex = (int) atomIndices.getParameterValue(p);
            return rewardRatesValues.getParameterValue((int) rewardRatesMapping.getParameterValue(stateIndex));
        } else {
            return super.getUntransformedBranchRate(tree, node);
        }
    }

    @Override
    public double getBranchRateDifferential(final Tree tree, final NodeRef node) {
        final int p = getParameterIndexFromNode(node);

        if (isOne(indicator.getParameterValue(p))) {
            return 0.0;
        } else {
            return super.getBranchRateDifferential(tree, node);
        }
    }

    @Override
    public double getBranchRateSecondDifferential(final Tree tree, final NodeRef node) {
        final int p = getParameterIndexFromNode(node);

        if (isOne(indicator.getParameterValue(p))) {
            return 0.0;
        } else {
            return super.getBranchRateSecondDifferential(tree, node);
        }
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        if (variable == indicator || variable == atomIndices ||
                variable == rewardRatesValues || variable == rewardRatesInternal ||
                variable == rewardRatesMapping) {
            fireModelChanged();
        } else {
            super.handleVariableChangedEvent(variable, index, type);
        }
    }

    private static boolean isOne(final double x) {
        final long r = Math.round(x);
        return Math.abs(x - r) <= 1e-9 && r == 1L;
    }
}
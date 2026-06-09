package dr.evomodel.treedatalikelihood.discrete;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evolution.tree.TreeTraitProvider;
import dr.evomodel.branchratemodel.RewardsAwareMixtureBranchRates;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.DataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.ProcessSimulation;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.hmc.NumericGradientStepSizeProvider;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.xml.Reportable;

import java.util.Arrays;

/**
 * Gradient of a dependent-process tree likelihood with respect to the free
 * reward-rate values used by atomic reward branches in a
 * {@link RewardsAwareMixtureBranchRates} model.
 *
 * Continuous branches are functions of total.rewards.cts and therefore have
 * zero derivative with respect to reward-rate values here. Atomic branches
 * contribute through the selected atom value.
 *
 * @author Filippo Monti
 */
public final class RewardsAwareMixtureBranchRatesRewardRateGradient
        implements GradientWrtParameterProvider, Reportable {

    private final TreeDataLikelihood treeDataLikelihood;
    private final Tree tree;
    private final BeagleDataLikelihoodDelegate likelihoodDelegate;
    private final RewardsAwareMixtureBranchRates branchRateModel;
    private final Parameter rewardRatesValues;
    private final Parameter rewardRatesVaryingValues;
    private final Parameter indicator;
    private final Parameter atomIndices;
    private final TreeTrait treeTraitProvider;
    private final Double tolerance;

    private double numericGradientStepSize =
            NumericGradientStepSizeProvider.StepSizeLevel.MEDIUM.getStepSizeRatio();

    public RewardsAwareMixtureBranchRatesRewardRateGradient(
            TreeDataLikelihood treeDataLikelihood,
            RewardsAwareMixtureBranchRates branchRateModel,
            Parameter rewardRatesValues,
            Parameter rewardRatesVaryingValues,
            Double tolerance) {

        if (treeDataLikelihood == null) {
            throw new IllegalArgumentException("treeDataLikelihood must be non-null");
        }
        if (branchRateModel == null) {
            throw new IllegalArgumentException("branchRateModel must be non-null");
        }
        if (rewardRatesValues == null) {
            throw new IllegalArgumentException("rewardRatesValues must be non-null");
        }
        if (rewardRatesVaryingValues == null) {
            throw new IllegalArgumentException("rewardRatesVaryingValues must be non-null");
        }
        if (rewardRatesVaryingValues.getDimension() != rewardRatesValues.getDimension() - 2) {
            throw new IllegalArgumentException(
                    "rewardRatesVaryingValues dimension must equal rewardRatesValues dimension - 2");
        }
        if (treeDataLikelihood.getBranchRateModel() != branchRateModel) {
            throw new IllegalArgumentException(
                    "TreeDataLikelihood must use the supplied RewardsAwareMixtureBranchRates model.");
        }

        final DataLikelihoodDelegate delegate = treeDataLikelihood.getDataLikelihoodDelegate();
        if (!(delegate instanceof BeagleDataLikelihoodDelegate)) {
            throw new IllegalArgumentException(
                    "RewardsAwareMixtureBranchRatesRewardRateGradient requires a BeagleDataLikelihoodDelegate");
        }

        this.treeDataLikelihood = treeDataLikelihood;
        this.tree = treeDataLikelihood.getTree();
        this.likelihoodDelegate = (BeagleDataLikelihoodDelegate) delegate;
        this.branchRateModel = branchRateModel;
        this.rewardRatesValues = rewardRatesValues;
        this.rewardRatesVaryingValues = rewardRatesVaryingValues;
        this.indicator = branchRateModel.getIndicator();
        this.atomIndices = branchRateModel.getAtomIndices();
        this.tolerance = tolerance;

        TreeTrait trait = treeDataLikelihood.getTreeTrait(DiscreteTraitBranchRateDelegate.getName(null));
        if (trait == null) {
            TreeTraitProvider traitProvider = new ProcessSimulation(
                    treeDataLikelihood,
                    new DiscreteTraitBranchRateDelegate(
                            "Sequence",
                            tree,
                            likelihoodDelegate,
                            branchRateModel));
            treeDataLikelihood.addTraits(traitProvider.getTreeTraits());
            trait = treeDataLikelihood.getTreeTrait(DiscreteTraitBranchRateDelegate.getName(null));
        }
        if (trait == null) {
            throw new IllegalStateException("Unable to create dependent branch-rate gradient trait.");
        }
        this.treeTraitProvider = trait;
    }

    @Override
    public Likelihood getLikelihood() {
        return treeDataLikelihood;
    }

    @Override
    public Parameter getParameter() {
        return rewardRatesVaryingValues;
    }

    @Override
    public int getDimension() {
        return rewardRatesVaryingValues.getDimension();
    }

    @Override
    public double[] getGradientLogDensity() {

        final double[] result = new double[getDimension()];
        final double[] branchGradient = (double[]) treeTraitProvider.getTrait(tree, null);

        int v = 0;
        for (int i = 0; i < tree.getNodeCount(); ++i) {
            final NodeRef node = tree.getNode(i);
            if (!tree.isRoot(node)) {
                final int parameterIndex = branchRateModel.getParameterIndexFromNode(node);
                if (isOne(indicator.getParameterValue(parameterIndex))) {
                    final int rewardRateIndex = checkedRewardRateIndex(parameterIndex);
                    final int varyingIndex = rewardRateIndex - 2;
                    if (varyingIndex >= 0) {
                        final double raw = rewardRatesValues.getParameterValue(rewardRateIndex);
                        result[varyingIndex] += branchGradient[v]
                                * tree.getBranchLength(node)
                                * branchRateModel.getTransform().differential(raw, tree, node);
                    }
                }
                v++;
            }
        }

        return result;
    }

    private int checkedRewardRateIndex(final int parameterIndex) {
        return mapAtomStateToRewardRateIndex(
                atomIndices,
                branchRateModel.getRewardRatesMapping(),
                rewardRatesValues,
                parameterIndex);
    }

    public static int mapAtomStateToRewardRateIndex(
            final Parameter atomIndices,
            final Parameter rewardRatesMapping,
            final Parameter rewardRatesValues,
            final int parameterIndex) {

        final double raw = atomIndices.getParameterValue(parameterIndex);
        final int stateIndex = (int) Math.round(raw);
        if (Math.abs(raw - stateIndex) > 1.0e-10 ||
                stateIndex < 0 ||
                stateIndex >= rewardRatesMapping.getDimension()) {
            throw new IllegalStateException(
                    "Invalid atom state index at branch parameter " + parameterIndex + ": " + raw);
        }

        final double mapped = rewardRatesMapping.getParameterValue(stateIndex);
        final int rewardRateIndex = (int) Math.round(mapped);
        if (Math.abs(mapped - rewardRateIndex) > 1.0e-10 ||
                rewardRateIndex < 0 ||
                rewardRateIndex >= rewardRatesValues.getDimension()) {
            throw new IllegalStateException(
                    "Invalid reward-rate mapping for atom state " + stateIndex + ": " + mapped);
        }
        return rewardRateIndex;
    }

    private static boolean isOne(final double x) {
        final long r = Math.round(x);
        return Math.abs(x - r) <= 1.0e-9 && r == 1L;
    }

    @Override
    public String getReport() {
        if (tolerance != null) {
            return "Dependent reward-rate gradient; check tolerance=" + tolerance + '\n' +
                    GradientWrtParameterProvider.getReportAndCheckForError(
                            this,
                            Double.NEGATIVE_INFINITY,
                            Double.POSITIVE_INFINITY,
                            tolerance);
        }

        return "Dependent reward-rate gradient (no check tolerance specified).\n" +
                "analytic: " + Arrays.toString(getGradientLogDensity()) + '\n';
    }

    @Override
    public double getNumericGradientStepSize() {
        return numericGradientStepSize;
    }

    @Override
    public void setNumericGradientStepSize(double ratio) {
        numericGradientStepSize = ratio;
    }
}

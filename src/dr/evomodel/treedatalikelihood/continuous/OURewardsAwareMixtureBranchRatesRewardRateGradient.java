package dr.evomodel.treedatalikelihood.continuous;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evomodel.branchratemodel.RewardsAwareMixtureBranchRates;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.preorder.BranchConditionalDistributionDelegate;
import dr.evomodel.treedatalikelihood.preorder.BranchSufficientStatistics;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.hmc.NumericGradientStepSizeProvider;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.xml.Reportable;

import java.util.Arrays;
import java.util.List;

import static dr.evomodel.treedatalikelihood.discrete.RewardsAwareMixtureBranchRatesRewardRateGradient.mapAtomStateToRewardRateIndex;

public final class OURewardsAwareMixtureBranchRatesRewardRateGradient
        implements GradientWrtParameterProvider, Reportable {

    private static final double ZERO_BRANCH_LENGTH_TOLERANCE = 1.0e-12;

    private final TreeDataLikelihood treeDataLikelihood;
    private final TreeTrait<List<BranchSufficientStatistics>> treeTraitProvider;
    private final Tree tree;
    private final RewardsAwareMixtureBranchRates branchRateModel;
    private final Parameter rewardRatesValues;
    private final Parameter rewardRatesVaryingValues;
    private final Parameter indicator;
    private final Parameter atomIndices;
    private final int nTraits;
    private final OUBranchTimeGradient branchTimeGradient;
    private final Double tolerance;

    private double numericGradientStepSize =
            NumericGradientStepSizeProvider.StepSizeLevel.MEDIUM.getStepSizeRatio();

    public OURewardsAwareMixtureBranchRatesRewardRateGradient(
            String traitName,
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
        if (!(treeDataLikelihood.getDataLikelihoodDelegate() instanceof ContinuousDataLikelihoodDelegate)) {
            throw new IllegalArgumentException(
                    "OURewardsAwareMixtureBranchRatesRewardRateGradient requires a ContinuousDataLikelihoodDelegate");
        }

        final ContinuousDataLikelihoodDelegate likelihoodDelegate =
                (ContinuousDataLikelihoodDelegate) treeDataLikelihood.getDataLikelihoodDelegate();

        this.treeDataLikelihood = treeDataLikelihood;
        this.tree = treeDataLikelihood.getTree();
        this.branchRateModel = branchRateModel;
        this.rewardRatesValues = rewardRatesValues;
        this.rewardRatesVaryingValues = rewardRatesVaryingValues;
        this.indicator = branchRateModel.getIndicator();
        this.atomIndices = branchRateModel.getAtomIndices();
        this.tolerance = tolerance;

        final String bcdName = BranchConditionalDistributionDelegate.getName(traitName);
        if (treeDataLikelihood.getTreeTrait(bcdName) == null) {
            likelihoodDelegate.addBranchConditionalDensityTrait(traitName);
        }

        @SuppressWarnings("unchecked")
        final TreeTrait<List<BranchSufficientStatistics>> unchecked =
                treeDataLikelihood.getTreeTrait(bcdName);
        this.treeTraitProvider = unchecked;
        if (treeTraitProvider == null) {
            throw new IllegalStateException("Unable to create branch conditional density trait: " + bcdName);
        }

        this.nTraits = likelihoodDelegate.getTraitCount();
        if (nTraits != 1) {
            throw new RuntimeException("Not yet implemented for >1 traits");
        }

        this.branchTimeGradient = new OUBranchTimeGradient(
                likelihoodDelegate.getTraitDim(),
                tree,
                likelihoodDelegate);
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

        treeDataLikelihood.makeDirty();

        final double[] result = new double[getDimension()];

        for (int i = 0; i < tree.getNodeCount(); ++i) {
            final NodeRef node = tree.getNode(i);
            if (tree.isRoot(node) || isZeroLengthBranch(node)) {
                continue;
            }

            final int parameterIndex = branchRateModel.getParameterIndexFromNode(node);
            if (!isOne(indicator.getParameterValue(parameterIndex))) {
                continue;
            }

            final int rewardRateIndex = mapAtomStateToRewardRateIndex(
                    atomIndices,
                    branchRateModel.getRewardRatesMapping(),
                    rewardRatesValues,
                    parameterIndex);
            final int varyingIndex = rewardRateIndex - 2;
            if (varyingIndex < 0) {
                continue;
            }

            final List<BranchSufficientStatistics> statisticsForNode = treeTraitProvider.getTrait(tree, node);
            if (statisticsForNode.size() != nTraits) {
                throw new IllegalStateException(
                        "Expected " + nTraits + " branch sufficient statistics but found " +
                                statisticsForNode.size());
            }

            final double raw = rewardRatesValues.getParameterValue(rewardRateIndex);
            final double rate = branchRateModel.getBranchRate(tree, node);
            if (rate == 0.0) {
                continue;
            }
            final double differential = branchRateModel.getTransform().differential(raw, tree, node);
            final double timeDerivative = getScaledTime(node) * differential / rate;

            for (int trait = 0; trait < nTraits; ++trait) {
                result[varyingIndex] +=
                        branchTimeGradient.getGradientForBranch(statisticsForNode.get(trait), node)[0] *
                                timeDerivative;
            }
        }

        return result;
    }

    private double getScaledTime(NodeRef node) {
        final ContinuousDataLikelihoodDelegate delegate =
                (ContinuousDataLikelihoodDelegate) treeDataLikelihood.getDataLikelihoodDelegate();
        final int matrixIndex = delegate.getDiffusionProcessDelegate()
                .getMatrixBufferOffsetIndex(node.getNumber());
        return delegate.getIntegrator().getBranchLength(matrixIndex);
    }

    private boolean isZeroLengthBranch(NodeRef node) {
        return Math.abs(tree.getBranchLength(node)) <= ZERO_BRANCH_LENGTH_TOLERANCE;
    }

    private static boolean isOne(final double x) {
        final long r = Math.round(x);
        return Math.abs(x - r) <= 1.0e-9 && r == 1L;
    }

    @Override
    public String getReport() {
        if (tolerance != null) {
            return "OU dependent reward-rate gradient; check tolerance=" + tolerance + '\n' +
                    GradientWrtParameterProvider.getReportAndCheckForError(
                            this,
                            Double.NEGATIVE_INFINITY,
                            Double.POSITIVE_INFINITY,
                            tolerance);
        }

        return "OU dependent reward-rate gradient (no check tolerance specified).\n" +
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

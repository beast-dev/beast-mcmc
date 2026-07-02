package dr.evomodel.treedatalikelihood.continuous;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evomodel.branchratemodel.ArbitraryBranchRates;
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

abstract class AbstractRewardsAwareMixtureBranchRatesLocationScaleGradient
        implements GradientWrtParameterProvider, Reportable {

    protected static final int LOCATION_INDEX = 0;
    protected static final int SCALE_INDEX = 1;

    private static final double ZERO_BRANCH_LENGTH_TOLERANCE = 1.0e-12;

    protected final TreeDataLikelihood treeDataLikelihood;
    protected final ContinuousDataLikelihoodDelegate likelihoodDelegate;
    protected final Tree tree;
    protected final RewardsAwareMixtureBranchRates branchRateModel;

    private final TreeTrait<List<BranchSufficientStatistics>> treeTraitProvider;
    private final Parameter locationScaleParameter;
    private final int nTraits;
    private final Double tolerance;
    private final String reportName;

    private double numericGradientStepSize =
            NumericGradientStepSizeProvider.StepSizeLevel.MEDIUM.getStepSizeRatio();

    AbstractRewardsAwareMixtureBranchRatesLocationScaleGradient(
            String reportName,
            String traitName,
            TreeDataLikelihood treeDataLikelihood,
            RewardsAwareMixtureBranchRates branchRateModel,
            Parameter locationScaleParameter,
            Double tolerance) {

        if (treeDataLikelihood == null) {
            throw new IllegalArgumentException("treeDataLikelihood must be non-null");
        }
        if (branchRateModel == null) {
            throw new IllegalArgumentException("branchRateModel must be non-null");
        }
        if (locationScaleParameter == null) {
            throw new IllegalArgumentException("locationScaleParameter must be non-null");
        }
        if (locationScaleParameter.getDimension() != 2) {
            throw new IllegalArgumentException("locationScaleParameter must have dimension 2");
        }
        if (treeDataLikelihood.getBranchRateModel() != branchRateModel) {
            throw new IllegalArgumentException(
                    "TreeDataLikelihood must use the supplied RewardsAwareMixtureBranchRates model.");
        }
        if (!(branchRateModel.getTransform() instanceof ArbitraryBranchRates.BranchRateTransform.LinearLocationScale)) {
            throw new IllegalArgumentException(
                    reportName + " requires linearLocationScale=\"true\".");
        }
        if (!(treeDataLikelihood.getDataLikelihoodDelegate() instanceof ContinuousDataLikelihoodDelegate)) {
            throw new IllegalArgumentException(
                    reportName + " requires a ContinuousDataLikelihoodDelegate");
        }

        this.reportName = reportName;
        this.treeDataLikelihood = treeDataLikelihood;
        this.likelihoodDelegate =
                (ContinuousDataLikelihoodDelegate) treeDataLikelihood.getDataLikelihoodDelegate();
        this.tree = treeDataLikelihood.getTree();
        this.branchRateModel = branchRateModel;
        this.locationScaleParameter = locationScaleParameter;
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
    }

    @Override
    public Likelihood getLikelihood() {
        return treeDataLikelihood;
    }

    @Override
    public Parameter getParameter() {
        return locationScaleParameter;
    }

    @Override
    public int getDimension() {
        return locationScaleParameter.getDimension();
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

            final double rate = branchRateModel.getBranchRate(tree, node);
            if (rate == 0.0) {
                continue;
            }

            final List<BranchSufficientStatistics> statisticsForNode = treeTraitProvider.getTrait(tree, node);
            if (statisticsForNode.size() != nTraits) {
                throw new IllegalStateException(
                        "Expected " + nTraits + " branch sufficient statistics but found " +
                                statisticsForNode.size());
            }

            final double rawReward = branchRateModel.getUntransformedBranchRate(tree, node);
            for (int trait = 0; trait < nTraits; ++trait) {
                accumulateBranchGradient(result, node, statisticsForNode.get(trait), rawReward, rate);
            }
        }

        return result;
    }

    protected abstract void accumulateBranchGradient(
            double[] result,
            NodeRef node,
            BranchSufficientStatistics statistics,
            double rawReward,
            double rate);

    private boolean isZeroLengthBranch(NodeRef node) {
        return Math.abs(tree.getBranchLength(node)) <= ZERO_BRANCH_LENGTH_TOLERANCE;
    }

    @Override
    public String getReport() {
        if (tolerance != null) {
            return reportName + "; check tolerance=" + tolerance + '\n' +
                    GradientWrtParameterProvider.getReportAndCheckForError(
                            this,
                            Double.NEGATIVE_INFINITY,
                            Double.POSITIVE_INFINITY,
                            tolerance);
        }

        return reportName + " (no check tolerance specified).\n" +
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

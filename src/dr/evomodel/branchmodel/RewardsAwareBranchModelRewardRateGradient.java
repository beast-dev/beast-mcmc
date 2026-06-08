package dr.evomodel.branchmodel;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.treedatalikelihood.DataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.DiscreteDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.hmc.NumericGradientStepSizeProvider;
import dr.inference.markovjumps.SericolaSeriesMarkovRewardFastModel;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.xml.Reportable;

import java.util.Arrays;

/**
 * Gradient of the log tree data likelihood with respect to the reward-rate
 * proportions used by a {@link RewardsAwareBranchModel}.
 *
 * For fixed mixture indicators and atom indices, atomic no-jump branches have
 * zero smooth derivative with respect to reward-rate values: their local factor
 * is exp(Q_jj t). Continuous branches contribute through the Sericola reward
 * density adjoint.
 *
 * Current assumptions match {@link RewardsAwareBranchModelGradient}: the
 * discrete delegate supplies first-pattern / first-category branch messages.
 *
 * @author Filippo Monti
 */
public final class RewardsAwareBranchModelRewardRateGradient
        implements GradientWrtParameterProvider, Reportable {

    private final TreeDataLikelihood treeDataLikelihood;
    private final Tree tree;
    private final DiscreteDataLikelihoodDelegate discreteDelegate;
    private final RewardsAwareBranchModel rewardsAwareBranchModel;
    private final BranchRateModel branchRateModel;
    private final SericolaSeriesMarkovRewardFastModel sericola;
    private final Parameter rewardRatesValues;
    private final Parameter indicator;
    private final int nstates;

    private final Double tolerance;
    private double numericGradientStepSize =
            NumericGradientStepSizeProvider.StepSizeLevel.MEDIUM.getStepSizeRatio();

    private final double[] prePartial;
    private final double[] prePartialAtNode;
    private final double[] postPartial;
    private final double[] densityAdjoint;
    private final double[] branchGradient;
    private final double[] gradientBuffer;

    public RewardsAwareBranchModelRewardRateGradient(
            TreeDataLikelihood treeDataLikelihood,
            RewardsAwareBranchModel rewardsAwareBranchModel,
            Parameter rewardRatesValues,
            Double tolerance) {

        if (treeDataLikelihood == null) {
            throw new IllegalArgumentException("treeDataLikelihood must be non-null");
        }
        if (rewardsAwareBranchModel == null) {
            throw new IllegalArgumentException("rewardsAwareBranchModel must be non-null");
        }
        if (rewardRatesValues == null) {
            throw new IllegalArgumentException("rewardRatesValues must be non-null");
        }

        this.treeDataLikelihood = treeDataLikelihood;
        this.rewardsAwareBranchModel = rewardsAwareBranchModel;
        this.tree = rewardsAwareBranchModel.getTree();
        this.branchRateModel = rewardsAwareBranchModel.getRateBranchModel();
        this.sericola = rewardsAwareBranchModel.getSericolaModel();
        this.rewardRatesValues = rewardRatesValues;
        this.indicator = rewardsAwareBranchModel.getIndicator();
        this.nstates = rewardsAwareBranchModel.getStateCount();
        this.tolerance = tolerance;

        if (treeDataLikelihood.getTree().getNodeCount() != tree.getNodeCount()) {
            throw new IllegalArgumentException(
                    "TreeDataLikelihood and RewardsAwareBranchModel must use trees with the same node count.");
        }

        final DataLikelihoodDelegate delegate = treeDataLikelihood.getDataLikelihoodDelegate();
        if (!(delegate instanceof DiscreteDataLikelihoodDelegate)) {
            throw new IllegalArgumentException(
                    "RewardsAwareBranchModelRewardRateGradient requires TreeDataLikelihood to use DiscreteDataLikelihoodDelegate");
        }
        this.discreteDelegate = (DiscreteDataLikelihoodDelegate) delegate;

        this.prePartial = new double[nstates];
        this.prePartialAtNode = new double[nstates];
        this.postPartial = new double[nstates];
        this.densityAdjoint = new double[nstates * nstates];
        this.branchGradient = new double[rewardRatesValues.getDimension()];
        this.gradientBuffer = new double[rewardRatesValues.getDimension()];
    }

    @Override
    public Likelihood getLikelihood() {
        return treeDataLikelihood;
    }

    @Override
    public Parameter getParameter() {
        return rewardRatesValues;
    }

    @Override
    public int getDimension() {
        return rewardRatesValues.getDimension();
    }

    @Override
    public double[] getGradientLogDensity() {

        discreteDelegate.updatePostOrdersFromTreeDataLikelihood(treeDataLikelihood);
        discreteDelegate.ensurePreOrderComputed();

        Arrays.fill(gradientBuffer, 0.0);

        for (int i = 0; i < tree.getNodeCount(); i++) {
            final NodeRef node = tree.getNode(i);
            if (tree.isRoot(node)) {
                continue;
            }

            final int nodeNumber = node.getNumber();
            final int parameterIndex = rewardsAwareBranchModel.getParameterIndexForNode(nodeNumber);
            if (indicator != null && isOne(indicator.getParameterValue(parameterIndex))) {
                continue;
            }

            accumulateContinuousBranchGradient(node, nodeNumber);
        }

        return gradientBuffer;
    }

    private void accumulateContinuousBranchGradient(final NodeRef node, final int nodeNumber) {

        discreteDelegate.getPreOrderAtBranchStartInto(nodeNumber, prePartial);
        discreteDelegate.getPreOrderAtBranchEndInto(nodeNumber, prePartialAtNode);
        discreteDelegate.getPostOrderAtBranchEndInto(nodeNumber, postPartial);

        final double branchLikelihood = dot(prePartialAtNode, postPartial);
        if (!(branchLikelihood > 0.0) ||
                Double.isNaN(branchLikelihood) ||
                Double.isInfinite(branchLikelihood)) {
            throw new IllegalStateException(
                    "Invalid branch likelihood normalization for node " + nodeNumber +
                            ": " + branchLikelihood +
                            ", prePartialAtNode=" + Arrays.toString(prePartialAtNode) +
                            ", postPartial=" + Arrays.toString(postPartial));
        }

        fillDensityAdjoint(prePartial, postPartial, 1.0 / branchLikelihood, densityAdjoint);

        final double branchLength = tree.getBranchLength(node);
        final double rewardProportion = branchRateModel.getBranchRate(tree, node);

        try {
            sericola.computePdfGradientWrtRewardRatesInto(
                    rewardProportion,
                    branchLength,
                    densityAdjoint,
                    branchGradient);
        } catch (RuntimeException e) {
            throw new RuntimeException(
                    "Failed to compute reward-rate gradient for node " + nodeNumber +
                            ", rewardProportion=" + rewardProportion +
                            ", branchLength=" + branchLength,
                    e);
        }

        for (int p = 0; p < gradientBuffer.length; p++) {
            gradientBuffer[p] += branchGradient[p];
        }
    }

    private void fillDensityAdjoint(
            final double[] pre,
            final double[] post,
            final double scale,
            final double[] out) {

        for (int i = 0; i < nstates; i++) {
            final double preI = pre[i] * scale;
            final int rowBase = i * nstates;
            for (int j = 0; j < nstates; j++) {
                out[rowBase + j] = preI * post[j];
            }
        }
    }

    private static double dot(final double[] x, final double[] y) {
        double acc = 0.0;
        for (int i = 0; i < x.length; i++) {
            acc += x[i] * y[i];
        }
        return acc;
    }

    private static boolean isOne(final double x) {
        final long r = Math.round(x);
        return Math.abs(x - r) <= 1.0e-9 && r == 1L;
    }

    @Override
    public String getReport() {
        if (tolerance != null) {
            return "Reward-aware branch-model reward-rate gradient; check tolerance=" + tolerance + '\n' +
                    GradientWrtParameterProvider.getReportAndCheckForError(
                            this,
                            Double.NEGATIVE_INFINITY,
                            Double.POSITIVE_INFINITY,
                            tolerance);
        }

        return "Reward-aware branch-model reward-rate gradient (no check tolerance specified).\n" +
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

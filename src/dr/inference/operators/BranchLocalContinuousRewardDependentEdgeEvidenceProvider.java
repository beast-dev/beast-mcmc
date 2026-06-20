package dr.inference.operators;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.RewardsAwareMixtureBranchRates;
import dr.evomodel.treedatalikelihood.DataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.continuous.DiffusionProcessDelegate;
import dr.evomodel.treedatalikelihood.continuous.cdi.ContinuousDiffusionIntegrator;
import dr.evomodel.treedatalikelihood.preorder.BranchConditionalDistributionDelegate;
import dr.evomodel.treedatalikelihood.preorder.BranchSufficientStatistics;
import dr.evomodel.treedatalikelihood.preorder.NormalSufficientStatistics;
import dr.evomodel.treedatalikelihood.preorder.PreOrderPrecision;
import dr.inference.model.Parameter;
import dr.math.matrixAlgebra.missingData.InversionResult;
import dr.math.matrixAlgebra.missingData.MissingOps;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

import java.util.List;

import static dr.evomodel.treedatalikelihood.continuous.cdi.ContinuousDiffusionIntegrator.LOG_SQRT_2_PI;

/**
 * Branch-local continuous-trait evidence for candidate reward values.
 *
 * This uses the existing continuous preorder/postorder branch conditional
 * statistics. The returned value may omit constants common to all candidate
 * rewards for the same branch, which is sufficient for Gibbs resampling.
 */
public final class BranchLocalContinuousRewardDependentEdgeEvidenceProvider
        implements RewardDependentEdgeEvidenceProvider {

    private static final String TRAIT_NAME = "rewardDependentEdgeEvidence";

    private final TreeDataLikelihood treeDataLikelihood;
    private final ContinuousDataLikelihoodDelegate likelihoodDelegate;
    private final TreeTrait<List<BranchSufficientStatistics>> treeTraitProvider;
    private final Tree tree;
    private final RewardsAwareMixtureBranchRates rewardBranchRates;
    private final DiffusionProcessDelegate diffusionProcessDelegate;
    private final ContinuousDiffusionIntegrator cdi;
    private final Parameter continuousRewards;
    private final int dim;
    private final int nTraits;

    private final int[] branchUpdateIndices = new int[1];
    private final double[] branchLengths = new double[1];

    private final Transition currentTransition;
    private final Transition candidateTransition;

    private final DenseMatrix64F parentMean;
    private final DenseMatrix64F parentVariance;
    private final DenseMatrix64F candidateMean;
    private final DenseMatrix64F candidateVariance;
    private final DenseMatrix64F totalVariance;
    private final DenseMatrix64F totalPrecision;
    private final DenseMatrix64F inverseActualization;
    private final DenseMatrix64F matrix0;
    private final DenseMatrix64F matrix1;
    private final DenseMatrix64F vector0;

    public BranchLocalContinuousRewardDependentEdgeEvidenceProvider(final TreeDataLikelihood treeDataLikelihood) {
        if (treeDataLikelihood == null) {
            throw new IllegalArgumentException("treeDataLikelihood must be non-null");
        }

        final DataLikelihoodDelegate delegate = treeDataLikelihood.getDataLikelihoodDelegate();
        if (!(delegate instanceof ContinuousDataLikelihoodDelegate)) {
            throw new IllegalArgumentException(
                    "Continuous reward evidence requires ContinuousDataLikelihoodDelegate, found " +
                            (delegate == null ? "null" : delegate.getClass().getName())
            );
        }

        final BranchRateModel branchRateModel = treeDataLikelihood.getBranchRateModel();
        if (!(branchRateModel instanceof RewardsAwareMixtureBranchRates)) {
            throw new IllegalArgumentException(
                    "Continuous reward evidence requires RewardsAwareMixtureBranchRates, found " +
                            branchRateModel.getClass().getName()
            );
        }

        this.treeDataLikelihood = treeDataLikelihood;
        this.likelihoodDelegate = (ContinuousDataLikelihoodDelegate) delegate;
        this.tree = treeDataLikelihood.getTree();
        this.rewardBranchRates = (RewardsAwareMixtureBranchRates) branchRateModel;
        this.diffusionProcessDelegate = likelihoodDelegate.getDiffusionProcessDelegate();
        this.cdi = likelihoodDelegate.getIntegrator();
        this.continuousRewards = rewardBranchRates.getRateParameter();
        this.dim = likelihoodDelegate.getTraitDim();
        this.nTraits = likelihoodDelegate.getTraitCount();

        if (nTraits != 1) {
            throw new UnsupportedOperationException(
                    "Branch-local continuous reward evidence currently supports exactly one trait.");
        }

        final String bcdName = BranchConditionalDistributionDelegate.getName(TRAIT_NAME);
        if (treeDataLikelihood.getTreeTrait(bcdName) == null) {
            likelihoodDelegate.addBranchConditionalDensityTrait(TRAIT_NAME);
        }

        @SuppressWarnings("unchecked")
        final TreeTrait<List<BranchSufficientStatistics>> unchecked =
                treeDataLikelihood.getTreeTrait(bcdName);
        this.treeTraitProvider = unchecked;
        if (treeTraitProvider == null) {
            throw new IllegalStateException("Unable to create branch conditional density trait: " + bcdName);
        }

        this.currentTransition = new Transition(dim);
        this.candidateTransition = new Transition(dim);
        this.parentMean = new DenseMatrix64F(dim, 1);
        this.parentVariance = new DenseMatrix64F(dim, dim);
        this.candidateMean = new DenseMatrix64F(dim, 1);
        this.candidateVariance = new DenseMatrix64F(dim, dim);
        this.totalVariance = new DenseMatrix64F(dim, dim);
        this.totalPrecision = new DenseMatrix64F(dim, dim);
        this.inverseActualization = new DenseMatrix64F(dim, dim);
        this.matrix0 = new DenseMatrix64F(dim, dim);
        this.matrix1 = new DenseMatrix64F(dim, dim);
        this.vector0 = new DenseMatrix64F(dim, 1);
    }

    @Override
    public void prepare() {
        treeDataLikelihood.makeDirty();
    }

    @Override
    public double logEvidence(final int branchNodeNumber, final double rawReward) {
        if (branchNodeNumber < 0 || branchNodeNumber >= tree.getNodeCount()) {
            throw new IllegalArgumentException("branchNodeNumber out of range: " + branchNodeNumber);
        }

        final NodeRef node = tree.getNode(branchNodeNumber);
        if (tree.isRoot(node)) {
            throw new IllegalArgumentException("Root node has no branch: " + branchNodeNumber);
        }

        final int parameterIndex = rewardBranchRates.getParameterIndexFromNode(node);
        if (parameterIndex < 0 || parameterIndex >= continuousRewards.getDimension()) {
            throw new IllegalArgumentException(
                    "Invalid continuous reward parameter index " + parameterIndex +
                            " for branch node " + branchNodeNumber);
        }

        final List<BranchSufficientStatistics> statisticsForNode = treeTraitProvider.getTrait(tree, node);
        if (statisticsForNode.size() != nTraits) {
            throw new IllegalStateException(
                    "Expected " + nTraits + " branch sufficient statistics but found " +
                            statisticsForNode.size());
        }

        final BranchSufficientStatistics statistics = statisticsForNode.get(0);
        requireDense(statistics.getBelow(), "below");
        requireDense(statistics.getAbove(), "above");

        final double currentLength = getCurrentBranchLength(node);
        try {
            readTransitionForLength(node, currentLength, currentTransition);
            recoverParentContext(statistics.getAbove(), currentTransition);

            readTransitionForLength(node, getCandidateBranchLength(node, rawReward), candidateTransition);
            propagateCandidateContext(candidateTransition);

            return logIntegratedProduct(statistics.getBelow());
        } finally {
            setBranchMatrix(node, currentLength);
        }
    }

    private void recoverParentContext(final NormalSufficientStatistics above,
                                      final Transition current) {
        CommonOps.subtract(above.getRawVariance(), current.variance, matrix0);

        CommonOps.invert(current.actualization, inverseActualization);
        CommonOps.mult(inverseActualization, matrix0, matrix1);
        CommonOps.multTransB(matrix1, inverseActualization, parentVariance);

        CommonOps.subtract(above.getRawMean(), current.displacement, vector0);
        CommonOps.mult(inverseActualization, vector0, parentMean);
    }

    private void propagateCandidateContext(final Transition candidate) {
        CommonOps.mult(candidate.actualization, parentMean, candidateMean);
        CommonOps.addEquals(candidateMean, candidate.displacement);

        CommonOps.mult(candidate.actualization, parentVariance, matrix0);
        CommonOps.multTransB(matrix0, candidate.actualization, candidateVariance);
        CommonOps.addEquals(candidateVariance, candidate.variance);
    }

    private double logIntegratedProduct(final NormalSufficientStatistics below) {
        CommonOps.add(candidateVariance, below.getRawVariance(), totalVariance);

        final InversionResult inversion = MissingOps.safeInvertVariance(totalVariance, totalPrecision, true);
        final int effectiveDim = inversion.getEffectiveDimension();
        if (effectiveDim == 0) {
            return 0.0;
        }

        double ss = 0.0;
        for (int row = 0; row < dim; row++) {
            final double rowDifference = candidateMean.unsafe_get(row, 0) - below.getMean(row);
            for (int col = 0; col < dim; col++) {
                final double colDifference = candidateMean.unsafe_get(col, 0) - below.getMean(col);
                ss += rowDifference * totalPrecision.unsafe_get(row, col) * colDifference;
            }
        }

        final double logEvidence =
                -effectiveDim * LOG_SQRT_2_PI -
                        0.5 * inversion.getLogDeterminant() -
                        0.5 * ss;
        return Double.isFinite(logEvidence) ? logEvidence : Double.NEGATIVE_INFINITY;
    }

    private void readTransitionForLength(final NodeRef node,
                                         final double branchLength,
                                         final Transition transition) {
        setBranchMatrix(node, branchLength);

        final int matrixIndex = diffusionProcessDelegate.getMatrixBufferOffsetIndex(node.getNumber());
        final int precisionIndex = diffusionProcessDelegate.getEigenBufferOffsetIndex(0);

        cdi.getBranchVariance(matrixIndex, precisionIndex, transition.varianceData);
        cdi.getBranchDisplacement(matrixIndex, transition.displacementData);
        readActualization(matrixIndex, transition);
    }

    private void setBranchMatrix(final NodeRef node, final double branchLength) {
        branchUpdateIndices[0] = node.getNumber();
        branchLengths[0] = branchLength;
        diffusionProcessDelegate.updateDiffusionMatrices(
                cdi,
                branchUpdateIndices,
                branchLengths,
                1,
                false
        );
    }

    private void readActualization(final int matrixIndex, final Transition transition) {
        CommonOps.setIdentity(transition.actualization);
        if (!diffusionProcessDelegate.hasActualization()) {
            return;
        }

        CommonOps.fill(transition.actualization, 0.0);
        if (diffusionProcessDelegate.hasDiagonalActualization()) {
            cdi.getBranchActualization(matrixIndex, transition.actualizationDiagonal);
            for (int i = 0; i < dim; i++) {
                transition.actualization.unsafe_set(i, i, transition.actualizationDiagonal[i]);
            }
        } else {
            cdi.getBranchActualization(matrixIndex, transition.actualizationData);
        }
    }

    private double getCurrentBranchLength(final NodeRef node) {
        return tree.getBranchLength(node) *
                rewardBranchRates.getBranchRate(tree, node) *
                likelihoodDelegate.getRateTransformationNormalization();
    }

    private double getCandidateBranchLength(final NodeRef node, final double rawReward) {
        return tree.getBranchLength(node) *
                rewardBranchRates.getBranchRateForRawReward(tree, node, rawReward) *
                likelihoodDelegate.getRateTransformationNormalization();
    }

    private static void requireDense(final NormalSufficientStatistics statistics,
                                     final String label) {
        if (!(statistics.getPrecision() instanceof PreOrderPrecision.Dense)) {
            throw new UnsupportedOperationException(
                    "Branch-local continuous reward evidence currently requires dense " +
                            label + " precision statistics.");
        }
    }

    private static final class Transition {
        private final double[] varianceData;
        private final double[] displacementData;
        private final double[] actualizationData;
        private final double[] actualizationDiagonal;

        private final DenseMatrix64F variance;
        private final DenseMatrix64F displacement;
        private final DenseMatrix64F actualization;

        private Transition(final int dim) {
            this.varianceData = new double[dim * dim];
            this.displacementData = new double[dim];
            this.actualizationData = new double[dim * dim];
            this.actualizationDiagonal = new double[dim];

            this.variance = DenseMatrix64F.wrap(dim, dim, varianceData);
            this.displacement = DenseMatrix64F.wrap(dim, 1, displacementData);
            this.actualization = DenseMatrix64F.wrap(dim, dim, actualizationData);
        }
    }
}

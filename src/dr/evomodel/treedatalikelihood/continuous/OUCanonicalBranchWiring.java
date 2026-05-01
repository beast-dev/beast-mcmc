package dr.evomodel.treedatalikelihood.continuous;

import dr.evomodel.treedatalikelihood.preorder.BranchSufficientStatistics;
import dr.evomodel.treedatalikelihood.preorder.MatrixSufficientStatistics;
import dr.evomodel.treedatalikelihood.preorder.NormalSufficientStatistics;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalBranchMessageContribution;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalLocalTransitionAdjoints;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalTransitionAdjointUtils;
import dr.evomodel.continuous.ou.OUProcessModel;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianTransition;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalNumerics;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalNumericsOptions;
import org.ejml.data.DenseMatrix64F;

/**
 * Local canonical OU branch factor assembly and adjoint extraction.
 */
public final class OUCanonicalBranchWiring {
    private static final CanonicalNumericsOptions NUMERICS_OPTIONS = CanonicalNumericsOptions.OU_TREE;

    private final OUGaussianBranchTransitionProvider branchTransitionProvider;
    private final OUProcessModel processModel;
    private final int dimension;
    private final CanonicalDebugOptions debugOptions;

    private final CanonicalGaussianTransition transition;
    private final CanonicalTransitionAdjointUtils.Workspace workspace;
    private final OUCanonicalTransitionState transitionState;
    private final OUCanonicalBranchContributionAssembler contributionAssembler;

    public OUCanonicalBranchWiring(final OUGaussianBranchTransitionProvider branchTransitionProvider) {
        this(branchTransitionProvider,
                CanonicalDebugOptions.fromSystemProperties(),
                CanonicalGradientFallbackPolicy.fromSystemProperties());
    }

    OUCanonicalBranchWiring(final OUGaussianBranchTransitionProvider branchTransitionProvider,
                            final CanonicalDebugOptions debugOptions,
                            final CanonicalGradientFallbackPolicy fallbackPolicy) {
        if (branchTransitionProvider == null) {
            throw new IllegalArgumentException("branchTransitionProvider must not be null");
        }
        if (debugOptions == null) {
            throw new IllegalArgumentException("debugOptions must not be null");
        }
        if (fallbackPolicy == null) {
            throw new IllegalArgumentException("fallbackPolicy must not be null");
        }
        this.branchTransitionProvider = branchTransitionProvider;
        this.processModel = branchTransitionProvider.getProcessModel();
        this.dimension = branchTransitionProvider.getStateDimension();
        this.debugOptions = debugOptions;
        this.transition = new CanonicalGaussianTransition(dimension);
        this.workspace = new CanonicalTransitionAdjointUtils.Workspace(dimension);
        this.transitionState = new OUCanonicalTransitionState(
                branchTransitionProvider,
                transition,
                NUMERICS_OPTIONS);
        final OUCanonicalParentAboveMessages parentAboveMessages =
                new OUCanonicalParentAboveMessages(transitionState);
        this.contributionAssembler = new OUCanonicalBranchContributionAssembler(
                dimension,
                fallbackPolicy,
                transition,
                transitionState,
                parentAboveMessages);
    }

    public OUGaussianBranchTransitionProvider getBranchTransitionProvider() {
        return branchTransitionProvider;
    }

    public OUProcessModel getProcessModel() {
        return processModel;
    }

    public int getDimension() {
        return dimension;
    }

    public void fillLocalAdjoints(final BranchSufficientStatistics statistics,
                                  final CanonicalLocalTransitionAdjoints out) {
        final CanonicalBranchMessageContribution contribution = prepareCanonicalContribution(statistics);
        CanonicalTransitionAdjointUtils.fillFromCanonicalTransition(
                transition, contribution, workspace, out);
    }

    public void fillLocalAdjoints(final double branchLength,
                                  final double[] optimum,
                                  final BranchSufficientStatistics statistics,
                                  final CanonicalLocalTransitionAdjoints out) {
        final CanonicalBranchMessageContribution contribution =
                prepareCanonicalContribution(branchLength, optimum, statistics);
        CanonicalTransitionAdjointUtils.fillFromCanonicalTransition(
                transition, contribution, workspace, out);
    }

    public CanonicalBranchMessageContribution prepareCanonicalContribution(final BranchSufficientStatistics statistics) {
        fillCanonicalTransition(statistics.getBranch());
        return prepareCanonicalContributionFromCurrentTransition(statistics);
    }

    public CanonicalBranchMessageContribution prepareCanonicalContribution(final double branchLength,
                                                                           final double[] optimum,
                                                                           final BranchSufficientStatistics statistics) {
        if (hasFiniteBranchTransitionStatistics(statistics.getBranch())) {
            fillCanonicalTransition(statistics.getBranch());
        } else {
            if (debugOptions.isNonFiniteBranchStatsEnabled()) {
                System.err.println("nonfiniteBranchStatsFallback branchLength=" + branchLength);
            }
            fillCanonicalTransitionFromKernel(branchLength, optimum);
        }
        return prepareCanonicalContributionFromCurrentTransition(statistics);
    }

    private CanonicalBranchMessageContribution prepareCanonicalContributionFromCurrentTransition(
            final BranchSufficientStatistics statistics) {
        return contributionAssembler.prepareFromCurrentTransition(statistics);
    }

    public double evaluateFrozenLocalLogFactor(final double branchLength,
                                               final double[] optimum,
                                               final NormalSufficientStatistics aboveParent,
                                               final NormalSufficientStatistics below) {
        return contributionAssembler.evaluateFrozenLocalLogFactor(
                branchLength, optimum, aboveParent, below);
    }

    public void fillTransitionMomentsFromKernel(final double branchLength,
                                                final double[] optimum,
                                                final double[] transitionMatrixOut,
                                                final double[] transitionOffsetOut,
                                                final double[] transitionCovarianceOut) {
        transitionState.fillTransitionMomentsFromKernel(
                branchLength,
                optimum,
                transitionMatrixOut,
                transitionOffsetOut,
                transitionCovarianceOut);
    }

    public void accumulateBranchMeanGradient(final MatrixSufficientStatistics branchStatistics,
                                             final double[] dLogL_df,
                                             final double[] out) {
        zero(out);
        final DenseMatrix64F actualization = branchStatistics.getRawActualization();
        for (int j = 0; j < dimension; ++j) {
            double sum = dLogL_df[j];
            for (int i = 0; i < dimension; ++i) {
                sum -= actualization.unsafe_get(i, j) * dLogL_df[i];
            }
            out[j] = sum;
        }
    }

    public double evaluateLocalSelectionScore(final double branchLength,
                                              final double[] optimum,
                                              final CanonicalLocalTransitionAdjoints adjoints) {
        final double[] transitionMatrix = new double[dimension * dimension];
        final double[] transitionCovariance = new double[dimension * dimension];
        final double[] transitionOffset = new double[dimension];
        fillTransitionMomentsFromKernel(branchLength, optimum, transitionMatrix, transitionOffset, transitionCovariance);
        if (!isFinite(transitionMatrix)) {
            throw new IllegalStateException("Non-finite transition matrix in local selection score");
        }
        if (!isFinite(transitionCovariance)) {
            throw new IllegalStateException("Non-finite transition covariance in local selection score");
        }
        if (!isFinite(transitionOffset)) {
            throw new IllegalStateException("Non-finite transition offset in local selection score");
        }

        double score = 0.0;
        double fContribution = 0.0;
        double matrixContribution = 0.0;
        double covarianceContribution = 0.0;
        for (int i = 0; i < dimension; ++i) {
            final double fTerm = adjoints.dLogL_df[i] * transitionOffset[i];
            score += fTerm;
            fContribution += fTerm;
            final int iOffset = i * dimension;
            for (int j = 0; j < dimension; ++j) {
                final double matrixTerm = adjoints.dLogL_dF[iOffset + j] * transitionMatrix[iOffset + j];
                final double covarianceTerm = adjoints.dLogL_dOmega[iOffset + j] * transitionCovariance[iOffset + j];
                score += matrixTerm;
                score += covarianceTerm;
                matrixContribution += matrixTerm;
                covarianceContribution += covarianceTerm;
            }
        }
        if (!Double.isFinite(score)) {
            throw new IllegalStateException(
                    "Non-finite accumulated local selection score"
                            + " fContribution=" + fContribution
                            + " matrixContribution=" + matrixContribution
                            + " covarianceContribution=" + covarianceContribution
                            + " max|dF|=" + maxAbs(adjoints.dLogL_dF)
                            + " max|dOmega|=" + maxAbs(adjoints.dLogL_dOmega)
                            + " max|df|=" + maxAbs(adjoints.dLogL_df)
                            + " max|F|=" + maxAbs(transitionMatrix)
                            + " max|Omega|=" + maxAbs(transitionCovariance)
                            + " max|f|=" + maxAbs(transitionOffset));
        }
        return score;
    }

    private void fillCanonicalTransition(final MatrixSufficientStatistics branchStatistics) {
        transitionState.fillFromStatistics(branchStatistics);
    }

    private void fillCanonicalTransitionFromKernel(final double branchLength,
                                                   final double[] optimum) {
        transitionState.fillFromKernel(branchLength, optimum);
    }

    private static void zero(final double[] vector) {
        for (int i = 0; i < vector.length; ++i) {
            vector[i] = 0.0;
        }
    }

    private static boolean isFinite(final double[] values) {
        return CanonicalNumerics.isFinite(values);
    }

    private static boolean isFinite(final DenseMatrix64F matrix) {
        return CanonicalNumerics.isFinite(matrix);
    }

    private static double maxAbs(final double[] values) {
        double max = 0.0;
        for (double value : values) {
            max = Math.max(max, Math.abs(value));
        }
        return max;
    }

    private static boolean hasFiniteBranchTransitionStatistics(final MatrixSufficientStatistics branch) {
        return isFinite(branch.getRawDisplacement())
                && isFinite(branch.getRawActualization())
                && isFinite(branch.getRawPrecision());
    }

}

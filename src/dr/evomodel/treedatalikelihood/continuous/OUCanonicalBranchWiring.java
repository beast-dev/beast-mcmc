package dr.evomodel.treedatalikelihood.continuous;

import dr.evomodel.treedatalikelihood.preorder.BranchSufficientStatistics;
import dr.evomodel.treedatalikelihood.preorder.MatrixSufficientStatistics;
import dr.evomodel.treedatalikelihood.preorder.NormalSufficientStatistics;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalBranchMessageContribution;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalLocalTransitionAdjoints;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalTransitionAdjointUtils;
import dr.evomodel.continuous.ou.OUProcessModel;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianState;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianTransition;
import dr.evomodel.treedatalikelihood.continuous.canonical.math.GaussianFormConverter;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalNumerics;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalNumericsOptions;
import dr.evomodel.treedatalikelihood.continuous.observationmodel.SufficientStatisticsTipObservationPattern;
import org.ejml.data.DenseMatrix64F;

/**
 * Self-contained local tree canonical wiring for one OU branch.
 *
 * <p>This class is intentionally narrow: given tree-side branch sufficient statistics,
 * it reconstructs the canonical branch factor that the time-series machinery expects,
 * handles the exact observed-tip conditioning limit, and produces local canonical
 * transition adjoints. That keeps the branch factor construction explicit and testable
 * outside the larger delegate/gradient classes.</p>
 */
public final class OUCanonicalBranchWiring {
    private static final CanonicalNumericsOptions NUMERICS_OPTIONS = CanonicalNumericsOptions.OU_TREE;

    private final OUGaussianBranchTransitionProvider branchTransitionProvider;
    private final OUProcessModel processModel;
    private final int dimension;
    private final CanonicalDebugOptions debugOptions;
    private final CanonicalGradientFallbackPolicy fallbackPolicy;

    private final CanonicalGaussianTransition transition;
    private final CanonicalGaussianState pairPosterior;
    private final CanonicalGaussianState currentPosterior;
    private final CanonicalBranchMessageContribution contribution;
    private final CanonicalTransitionAdjointUtils.Workspace workspace;
    private final OUCanonicalTransitionState transitionState;
    private final OUCanonicalFrozenLocalFactorEvaluator frozenLocalFactorEvaluator;

    private final double[] aboveInformation;
    private final double[] belowInformation;
    private final SufficientStatisticsTipObservationPattern observationPattern;
    private final double[][] reducedPrecisionScratch;
    private final double[][] reducedCovarianceScratch;
    private final double[] reducedInformationScratch;
    private final double[] reducedMeanScratch;
    private final double[] currentMeanScratch;
    private final double[][] currentCovarianceScratch;
    private final double[][] reducedCholeskyScratch;
    private final double[][] reducedLowerInverseScratch;

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
        this.fallbackPolicy = fallbackPolicy;
        this.transition = new CanonicalGaussianTransition(dimension);
        this.pairPosterior = new CanonicalGaussianState(2 * dimension);
        this.currentPosterior = new CanonicalGaussianState(dimension);
        this.contribution = new CanonicalBranchMessageContribution(dimension);
        this.workspace = new CanonicalTransitionAdjointUtils.Workspace(dimension);
        this.observationPattern = new SufficientStatisticsTipObservationPattern(dimension);
        this.transitionState = new OUCanonicalTransitionState(
                branchTransitionProvider,
                transition,
                NUMERICS_OPTIONS);
        this.frozenLocalFactorEvaluator = new OUCanonicalFrozenLocalFactorEvaluator(
                dimension,
                transition,
                observationPattern);
        this.aboveInformation = new double[dimension];
        this.belowInformation = new double[dimension];
        final int maxReducedDimension = 2 * dimension;
        this.reducedPrecisionScratch = new double[maxReducedDimension][maxReducedDimension];
        this.reducedCovarianceScratch = new double[maxReducedDimension][maxReducedDimension];
        this.reducedInformationScratch = new double[maxReducedDimension];
        this.reducedMeanScratch = new double[maxReducedDimension];
        this.currentMeanScratch = new double[dimension];
        this.currentCovarianceScratch = new double[dimension][dimension];
        this.reducedCholeskyScratch = new double[maxReducedDimension][maxReducedDimension];
        this.reducedLowerInverseScratch = new double[maxReducedDimension][maxReducedDimension];
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
        prepareCanonicalContribution(statistics);
        CanonicalTransitionAdjointUtils.fillFromCanonicalTransition(
                transition, contribution, workspace, out);
    }

    public void fillLocalAdjoints(final double branchLength,
                                  final double[] optimum,
                                  final BranchSufficientStatistics statistics,
                                  final CanonicalLocalTransitionAdjoints out) {
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
        final int observedCount = fallbackPolicy.isExactTipShortcutDisabled()
                ? 0
                : classifyExactObservationPattern(statistics.getBelow());
        if (observedCount > 0) {
            fillContributionForPartiallyObservedTip(
                    statistics.getAbove(),
                    statistics.getAboveParent(),
                    statistics.getBelow(),
                    observedCount);
        } else {
            fillPairPosterior(
                    statistics.getAbove(),
                    statistics.getAboveParent(),
                    statistics.getBelow());
            fillContributionFromPairPosterior();
        }
        return contribution;
    }

    public double evaluateFrozenLocalLogFactor(final double branchLength,
                                               final double[] optimum,
                                               final NormalSufficientStatistics aboveChild,
                                               final NormalSufficientStatistics below) {
        transitionState.fillFromKernel(branchLength, optimum);
        final DenseMatrix64F parentAbovePrecision = recoverParentAbovePrecision(aboveChild);
        final DenseMatrix64F parentAboveMean = transitionState.getAboveParentMeanVector();
        return frozenLocalFactorEvaluator.evaluate(parentAbovePrecision, parentAboveMean, below);
    }

    public double evaluateFrozenLocalLogFactor(final double branchLength,
                                               final double[] optimum,
                                               final NormalSufficientStatistics aboveChild,
                                               final NormalSufficientStatistics aboveParent,
                                               final NormalSufficientStatistics below) {
        transitionState.fillFromKernel(branchLength, optimum);
        final DenseMatrix64F parentAbovePrecision = recoverOrUseParentAbovePrecision(aboveChild, aboveParent);
        final DenseMatrix64F parentAboveMean = transitionState.getAboveParentMeanVector();
        return frozenLocalFactorEvaluator.evaluate(parentAbovePrecision, parentAboveMean, below);
    }

    public void fillTransitionMomentsFromKernel(final double branchLength,
                                                final double[] optimum,
                                                final double[][] transitionMatrixOut,
                                                final double[] transitionOffsetOut,
                                                final double[][] transitionCovarianceOut) {
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
        final double[][] transitionMatrix = new double[dimension][dimension];
        final double[][] transitionCovariance = new double[dimension][dimension];
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
                final double matrixTerm = adjoints.dLogL_dF[iOffset + j] * transitionMatrix[i][j];
                final double covarianceTerm = adjoints.dLogL_dOmega[iOffset + j] * transitionCovariance[i][j];
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

    private static boolean isFinite(final double[][] values) {
        return CanonicalNumerics.isFinite(values);
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

    private static double maxAbs(final double[][] values) {
        double max = 0.0;
        for (double[] row : values) {
            max = Math.max(max, maxAbs(row));
        }
        return max;
    }

    private void fillPairPosterior(final NormalSufficientStatistics above,
                                   final NormalSufficientStatistics aboveParent,
                                   final NormalSufficientStatistics below) {
        final DenseMatrix64F abovePrecision = recoverOrUseParentAbovePrecision(above, aboveParent);
        final DenseMatrix64F aboveMean = transitionState.getAboveParentMeanVector();
        final DenseMatrix64F belowPrecision = below.getRawPrecision();
        final DenseMatrix64F belowMean = below.getRawMean();

        fillInformation(abovePrecision, aboveMean, aboveInformation);
        fillInformation(belowPrecision, belowMean, belowInformation);

        final int pairDimension = pairPosterior.getDimension();
        for (int i = 0; i < dimension; ++i) {
            pairPosterior.information[i] = transition.informationX[i] + aboveInformation[i];
            pairPosterior.information[dimension + i] = transition.informationY[i] + belowInformation[i];
            final int pairRowOffset = i * pairDimension;
            final int transitionRowOffset = i * dimension;
            final int lowerRowOffset = (dimension + i) * pairDimension;
            for (int j = 0; j < dimension; ++j) {
                pairPosterior.precision[pairRowOffset + j] =
                        transition.precisionXX[transitionRowOffset + j] + abovePrecision.unsafe_get(i, j);
                pairPosterior.precision[pairRowOffset + dimension + j] =
                        transition.precisionXY[transitionRowOffset + j];
                pairPosterior.precision[lowerRowOffset + j] =
                        transition.precisionYX[transitionRowOffset + j];
                pairPosterior.precision[lowerRowOffset + dimension + j] =
                        transition.precisionYY[transitionRowOffset + j] + belowPrecision.unsafe_get(i, j);
            }
        }
        pairPosterior.logNormalizer = 0.0;

        if (!isFinite(pairPosterior.precision) || !isFinite(pairPosterior.information)) {
            throw new IllegalStateException(
                    "Non-finite pair posterior before canonical->moments conversion"
                            + "; abovePrecisionFinite=" + isFinite(abovePrecision)
                            + "; belowPrecisionFinite=" + isFinite(belowPrecision)
                            + "; belowMeanFinite=" + isFinite(belowMean)
                            + "; transitionXXFinite=" + isFinite(transition.precisionXX)
                            + "; transitionXYFinite=" + isFinite(transition.precisionXY)
                            + "; transitionYYFinite=" + isFinite(transition.precisionYY)
                            + "; transitionInfoXFinite=" + isFinite(transition.informationX)
                            + "; transitionInfoYFinite=" + isFinite(transition.informationY));
        }
    }

    private void fillContributionForObservedTip(final NormalSufficientStatistics above,
                                                final NormalSufficientStatistics aboveParent,
                                                final NormalSufficientStatistics below) {
        final DenseMatrix64F abovePrecision = recoverOrUseParentAbovePrecision(above, aboveParent);
        final DenseMatrix64F aboveMean = transitionState.getAboveParentMeanVector();
        final DenseMatrix64F observedChild = below.getRawMean();

        fillInformation(abovePrecision, aboveMean, aboveInformation);

        for (int i = 0; i < dimension; ++i) {
            double info = transition.informationX[i] + aboveInformation[i];
            final int iOffset = i * dimension;
            for (int j = 0; j < dimension; ++j) {
                info -= transition.precisionXY[iOffset + j] * observedChild.unsafe_get(j, 0);
                currentPosterior.precision[iOffset + j] =
                        transition.precisionXX[iOffset + j] + abovePrecision.unsafe_get(i, j);
            }
            currentPosterior.information[i] = info;
        }
        currentPosterior.logNormalizer = 0.0;

        GaussianFormConverter.fillMomentsFromState(currentPosterior, currentMeanScratch, currentCovarianceScratch);

        OUCanonicalBranchContributionBuilder.fillObservedTipContribution(
                dimension,
                currentMeanScratch,
                currentCovarianceScratch,
                observedChild,
                contribution);
    }

    private void fillContributionForPartiallyObservedTip(final NormalSufficientStatistics above,
                                                         final NormalSufficientStatistics aboveParent,
                                                         final NormalSufficientStatistics below,
                                                         final int observedCount) {
        final DenseMatrix64F abovePrecision = recoverOrUseParentAbovePrecision(above, aboveParent);
        final DenseMatrix64F aboveMean = transitionState.getAboveParentMeanVector();
        final DenseMatrix64F observedChild = below.getRawMean();

        fillInformation(abovePrecision, aboveMean, aboveInformation);
        final int missingCount = collectObservationPartition(observedCount);
        final int reducedDimension = dimension + missingCount;
        final double[][] reducedPrecision = reducedPrecisionScratch;
        final double[][] reducedCovariance = reducedCovarianceScratch;
        final double[] reducedInformation = reducedInformationScratch;
        final double[] reducedMean = reducedMeanScratch;

        fillReducedCanonicalSystem(abovePrecision, observedChild, missingCount, reducedPrecision, reducedInformation);
        invertSymmetricPositiveDefinite(reducedPrecision, reducedDimension, reducedCovariance);
        for (int i = 0; i < reducedDimension; ++i) {
            double sum = 0.0;
            for (int j = 0; j < reducedDimension; ++j) {
                sum += reducedCovariance[i][j] * reducedInformation[j];
            }
            reducedMean[i] = sum;
        }

        OUCanonicalBranchContributionBuilder.fillPartiallyObservedTipContribution(
                dimension,
                observationPattern,
                observedChild,
                reducedMean,
                reducedCovariance,
                contribution);
    }

    private void fillContributionFromPairPosterior() {
        fillMomentsFromCanonicalLocally(pairPosterior, reducedMeanScratch, reducedCovarianceScratch);

        OUCanonicalBranchContributionBuilder.fillPairPosteriorContribution(
                dimension,
                reducedMeanScratch,
                reducedCovarianceScratch,
                contribution);
    }

    private DenseMatrix64F recoverOrUseParentAbovePrecision(final NormalSufficientStatistics aboveChild,
                                                            final NormalSufficientStatistics aboveParent) {
        return transitionState.recoverOrUseParentAbovePrecision(aboveChild, aboveParent);
    }

    private DenseMatrix64F recoverParentAbovePrecision(final NormalSufficientStatistics above) {
        return transitionState.recoverParentAbovePrecision(above);
    }

    private static boolean isFinite(final DenseMatrix64F matrix) {
        return CanonicalNumerics.isFinite(matrix);
    }

    private static void fillInformation(final DenseMatrix64F precision,
                                        final DenseMatrix64F mean,
                                        final double[] out) {
        final int d = out.length;
        for (int i = 0; i < d; ++i) {
            double sum = 0.0;
            for (int j = 0; j < d; ++j) {
                sum += precision.unsafe_get(i, j) * mean.unsafe_get(j, 0);
            }
            out[i] = sum;
        }
    }

    private int classifyExactObservationPattern(final NormalSufficientStatistics below) {
        return observationPattern.classify(below);
    }

    private int collectObservationPartition(final int observedCount) {
        return observationPattern.collectPartition(observedCount);
    }

    private void fillReducedCanonicalSystem(final DenseMatrix64F abovePrecision,
                                            final DenseMatrix64F observedChild,
                                            final int missingCount,
                                            final double[][] reducedPrecision,
                                            final double[] reducedInformation) {
        for (int i = 0; i < dimension; ++i) {
            double information = transition.informationX[i] + aboveInformation[i];
            final int iOffset = i * dimension;
            for (int j = 0; j < dimension; ++j) {
                reducedPrecision[i][j] = transition.precisionXX[iOffset + j] + abovePrecision.unsafe_get(i, j);
            }
            for (int observed = 0; observed < dimension - missingCount; ++observed) {
                final int observedIndex = observationPattern.observedIndex(observed);
                information -= transition.precisionXY[iOffset + observedIndex] * observedChild.unsafe_get(observedIndex, 0);
            }
            reducedInformation[i] = information;
            for (int missing = 0; missing < missingCount; ++missing) {
                reducedPrecision[i][dimension + missing] =
                        transition.precisionXY[iOffset + observationPattern.missingIndex(missing)];
            }
        }

        for (int missing = 0; missing < missingCount; ++missing) {
            final int childIndex = observationPattern.missingIndex(missing);
            final int row = dimension + missing;
            double information = transition.informationY[childIndex];
            final int childOffset = childIndex * dimension;
            for (int observed = 0; observed < dimension - missingCount; ++observed) {
                final int observedIndex = observationPattern.observedIndex(observed);
                information -= transition.precisionYY[childOffset + observedIndex] * observedChild.unsafe_get(observedIndex, 0);
            }
            reducedInformation[row] = information;
            for (int j = 0; j < dimension; ++j) {
                reducedPrecision[row][j] = transition.precisionYX[childOffset + j];
            }
            for (int otherMissing = 0; otherMissing < missingCount; ++otherMissing) {
                reducedPrecision[row][dimension + otherMissing] =
                        transition.precisionYY[childOffset + observationPattern.missingIndex(otherMissing)];
            }
        }
    }

    private void fillMomentsFromCanonicalLocally(final CanonicalGaussianState state,
                                                 final double[] meanOut,
                                                 final double[][] covarianceOut) {
        final int dimensionUsed = state.getDimension();
        invertSymmetricPositiveDefinite(state.precision, dimensionUsed, covarianceOut);
        for (int i = 0; i < dimensionUsed; ++i) {
            double sum = 0.0;
            for (int j = 0; j < dimensionUsed; ++j) {
                sum += covarianceOut[i][j] * state.information[j];
            }
            meanOut[i] = sum;
        }
    }

    private double invertSymmetricPositiveDefinite(final double[][] matrix,
                                                   final int dimensionUsed,
                                                   final double[][] inverseOut) {
        return CanonicalNumerics.invertSymmetricPositiveDefinite(
                matrix,
                dimensionUsed,
                inverseOut,
                reducedCholeskyScratch,
                reducedPrecisionScratch,
                reducedLowerInverseScratch);
    }

    private double invertSymmetricPositiveDefinite(final double[] matrix,
                                                   final int dimensionUsed,
                                                   final double[][] inverseOut) {
        final double[][] square = reducedPrecisionScratch;
        for (int i = 0; i < dimensionUsed; ++i) {
            System.arraycopy(matrix, i * dimensionUsed, square[i], 0, dimensionUsed);
        }
        return invertSymmetricPositiveDefinite(square, dimensionUsed, inverseOut);
    }

}

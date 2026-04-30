package dr.evomodel.treedatalikelihood.continuous;

import dr.evomodel.treedatalikelihood.preorder.BranchSufficientStatistics;
import dr.evomodel.treedatalikelihood.preorder.MatrixSufficientStatistics;
import dr.evomodel.treedatalikelihood.preorder.NormalSufficientStatistics;
import dr.evomodel.treedatalikelihood.continuous.gaussian.message.CanonicalBranchMessageContribution;
import dr.evomodel.treedatalikelihood.continuous.gaussian.message.CanonicalLocalTransitionAdjoints;
import dr.evomodel.treedatalikelihood.continuous.gaussian.message.CanonicalTransitionAdjointUtils;
import dr.evomodel.continuous.ou.OUProcessModel;
import dr.evomodel.treedatalikelihood.continuous.gaussian.CanonicalGaussianState;
import dr.evomodel.treedatalikelihood.continuous.gaussian.CanonicalGaussianTransition;
import dr.evomodel.treedatalikelihood.continuous.gaussian.CanonicalGaussianUtils;
import dr.evomodel.treedatalikelihood.continuous.gaussian.message.CanonicalNumerics;
import dr.evomodel.treedatalikelihood.continuous.gaussian.message.CanonicalNumericsOptions;
import org.ejml.data.DenseMatrix64F;

/**
 * Self-contained local tree-to-time-series canonical wiring for one OU branch.
 *
 * <p>This class is intentionally narrow: given tree-side branch sufficient statistics,
 * it reconstructs the canonical branch factor that the time-series machinery expects,
 * handles the exact observed-tip conditioning limit, and produces local canonical
 * transition adjoints. That keeps the branch factor construction explicit and testable
 * outside the larger delegate/gradient classes.</p>
 */
public final class OUCanonicalBranchWiring {
    private static final String DISABLE_EXACT_TIP_SHORTCUT_PROPERTY =
            "beast.experimental.disableExactTipShortcut";
    private static final String NONFINITE_BRANCH_STATS_DEBUG_PROPERTY =
            "beast.debug.ou.nonfiniteBranchStats";

    private static final CanonicalNumericsOptions NUMERICS_OPTIONS = CanonicalNumericsOptions.OU_TREE;

    private final TimeSeriesOUGaussianBranchTransitionProvider branchTransitionProvider;
    private final OUProcessModel processModel;
    private final int dimension;

    private final CanonicalGaussianTransition transition;
    private final CanonicalGaussianState pairPosterior;
    private final CanonicalGaussianState currentPosterior;
    private final CanonicalBranchMessageContribution contribution;
    private final CanonicalTransitionAdjointUtils.Workspace workspace;

    private final double[] aboveInformation;
    private final double[] belowInformation;
    private final boolean[] exactObservationMask;
    private final int[] observedIndexScratch;
    private final int[] missingIndexScratch;
    private final int[] reducedIndexByTrait;
    private final DenseMatrix64F actualizationMatrix;
    private final DenseMatrix64F displacementVector;
    private final DenseMatrix64F branchPrecisionMatrix;
    private final DenseMatrix64F branchCovarianceMatrix;
    private final OUCanonicalParentAboveResolver parentAboveResolver;
    private final DenseMatrix64F scratchMatrix;
    private final DenseMatrix64F secondaryScratchMatrix;
    private final double[][] reducedPrecisionScratch;
    private final double[][] reducedCovarianceScratch;
    private final double[] reducedInformationScratch;
    private final double[] reducedMeanScratch;
    private final double[] currentMeanScratch;
    private final double[][] currentCovarianceScratch;
    private final double[][] reducedCholeskyScratch;
    private final double[][] reducedLowerInverseScratch;
    private final double[][] transitionMatrixScratch;
    private final double[][] transitionCovarianceScratch;
    private final double[] transitionOffsetScratch;

    public OUCanonicalBranchWiring(final TimeSeriesOUGaussianBranchTransitionProvider branchTransitionProvider) {
        if (branchTransitionProvider == null) {
            throw new IllegalArgumentException("branchTransitionProvider must not be null");
        }
        this.branchTransitionProvider = branchTransitionProvider;
        this.processModel = branchTransitionProvider.getProcessModel();
        this.dimension = branchTransitionProvider.getStateDimension();
        this.transition = new CanonicalGaussianTransition(dimension);
        this.pairPosterior = new CanonicalGaussianState(2 * dimension);
        this.currentPosterior = new CanonicalGaussianState(dimension);
        this.contribution = new CanonicalBranchMessageContribution(dimension);
        this.workspace = new CanonicalTransitionAdjointUtils.Workspace(dimension);
        this.aboveInformation = new double[dimension];
        this.belowInformation = new double[dimension];
        this.exactObservationMask = new boolean[dimension];
        this.observedIndexScratch = new int[dimension];
        this.missingIndexScratch = new int[dimension];
        this.reducedIndexByTrait = new int[dimension];
        this.actualizationMatrix = new DenseMatrix64F(dimension, dimension);
        this.displacementVector = new DenseMatrix64F(dimension, 1);
        this.branchPrecisionMatrix = new DenseMatrix64F(dimension, dimension);
        this.branchCovarianceMatrix = new DenseMatrix64F(dimension, dimension);
        this.parentAboveResolver = new OUCanonicalParentAboveResolver(
                dimension,
                actualizationMatrix,
                displacementVector,
                branchPrecisionMatrix,
                branchCovarianceMatrix,
                NUMERICS_OPTIONS);
        this.scratchMatrix = new DenseMatrix64F(dimension, dimension);
        this.secondaryScratchMatrix = new DenseMatrix64F(dimension, dimension);
        final int maxReducedDimension = 2 * dimension;
        this.reducedPrecisionScratch = new double[maxReducedDimension][maxReducedDimension];
        this.reducedCovarianceScratch = new double[maxReducedDimension][maxReducedDimension];
        this.reducedInformationScratch = new double[maxReducedDimension];
        this.reducedMeanScratch = new double[maxReducedDimension];
        this.currentMeanScratch = new double[dimension];
        this.currentCovarianceScratch = new double[dimension][dimension];
        this.reducedCholeskyScratch = new double[maxReducedDimension][maxReducedDimension];
        this.reducedLowerInverseScratch = new double[maxReducedDimension][maxReducedDimension];
        this.transitionMatrixScratch = new double[dimension][dimension];
        this.transitionCovarianceScratch = new double[dimension][dimension];
        this.transitionOffsetScratch = new double[dimension];
    }

    public TimeSeriesOUGaussianBranchTransitionProvider getBranchTransitionProvider() {
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
            if (Boolean.getBoolean(NONFINITE_BRANCH_STATS_DEBUG_PROPERTY)) {
                System.err.println("nonfiniteBranchStatsFallback branchLength=" + branchLength);
            }
            fillCanonicalTransitionFromKernel(branchLength, optimum);
        }
        return prepareCanonicalContributionFromCurrentTransition(statistics);
    }

    private CanonicalBranchMessageContribution prepareCanonicalContributionFromCurrentTransition(
            final BranchSufficientStatistics statistics) {
        final boolean disableExactTipShortcut = Boolean.getBoolean(DISABLE_EXACT_TIP_SHORTCUT_PROPERTY);
        final int observedCount = disableExactTipShortcut
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
        fillCanonicalTransitionFromKernel(branchLength, optimum);
        final DenseMatrix64F parentAbovePrecision = recoverParentAbovePrecision(aboveChild);
        final DenseMatrix64F parentAboveMean = parentAboveResolver.getAboveParentMeanVector();
        return evaluateFrozenLocalLogFactorWithResolvedParentPrecision(parentAbovePrecision, parentAboveMean, below);
    }

    public double evaluateFrozenLocalLogFactor(final double branchLength,
                                               final double[] optimum,
                                               final NormalSufficientStatistics aboveChild,
                                               final NormalSufficientStatistics aboveParent,
                                               final NormalSufficientStatistics below) {
        fillCanonicalTransitionFromKernel(branchLength, optimum);
        final DenseMatrix64F parentAbovePrecision = recoverOrUseParentAbovePrecision(aboveChild, aboveParent);
        final DenseMatrix64F parentAboveMean = parentAboveResolver.getAboveParentMeanVector();
        return evaluateFrozenLocalLogFactorWithResolvedParentPrecision(parentAbovePrecision, parentAboveMean, below);
    }

    private double evaluateFrozenLocalLogFactorWithResolvedParentPrecision(
            final DenseMatrix64F parentAbovePrecision,
            final DenseMatrix64F parentAboveMean,
            final NormalSufficientStatistics below) {
        fillInformation(parentAbovePrecision, parentAboveMean, aboveInformation);

        final int observedCount = classifyExactObservationPattern(below);
        if (observedCount == dimension) {
            final DenseMatrix64F observed = below.getRawMean();
            final double[][] precision = new double[dimension][dimension];
            final double[] information = new double[dimension];
            for (int i = 0; i < dimension; ++i) {
                double info = transition.informationX[i] + aboveInformation[i];
                final int iOffset = i * dimension;
                for (int j = 0; j < dimension; ++j) {
                    precision[i][j] = transition.precisionXX[iOffset + j] + parentAbovePrecision.unsafe_get(i, j);
                    info -= transition.precisionXY[iOffset + j] * observed.unsafe_get(j, 0);
                }
                information[i] = info;
            }
            return normalizedLogNormalizer(precision, information, dimension)
                    - transition.logNormalizer
                    - 0.5 * quadraticForm(transition.precisionYY, observed)
                    + dot(transition.informationY, observed);
        }

        if (observedCount > 0) {
            return evaluateFrozenPartiallyObservedLocalLogFactor(
                    parentAbovePrecision,
                    below.getRawMean(),
                    observedCount);
        }

        fillInformation(below.getRawPrecision(), below.getRawMean(), belowInformation);
        final int doubled = 2 * dimension;
        final double[][] precision = new double[doubled][doubled];
        final double[] information = new double[doubled];
        for (int i = 0; i < dimension; ++i) {
            information[i] = transition.informationX[i] + aboveInformation[i];
            information[dimension + i] = transition.informationY[i] + belowInformation[i];
            final int iOffset = i * dimension;
            for (int j = 0; j < dimension; ++j) {
                precision[i][j] = transition.precisionXX[iOffset + j] + parentAbovePrecision.unsafe_get(i, j);
                precision[i][dimension + j] = transition.precisionXY[iOffset + j];
                precision[dimension + i][j] = transition.precisionYX[iOffset + j];
                precision[dimension + i][dimension + j] =
                        transition.precisionYY[iOffset + j] + below.getRawPrecision().unsafe_get(i, j);
            }
        }
        return normalizedLogNormalizer(precision, information, doubled) - transition.logNormalizer;
    }

    public void fillTransitionMomentsFromKernel(final double branchLength,
                                                final double[] optimum,
                                                final double[][] transitionMatrixOut,
                                                final double[] transitionOffsetOut,
                                                final double[][] transitionCovarianceOut) {
        branchTransitionProvider.fillBranchTransitionMatrix(branchLength, transitionMatrixOut);
        branchTransitionProvider.fillBranchTransitionCovariance(branchLength, transitionCovarianceOut);
        fillTransitionOffset(transitionMatrixOut, optimum, transitionOffsetOut);
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
        final DenseMatrix64F actualization = branchStatistics.getRawActualization();
        final DenseMatrix64F displacement = branchStatistics.getRawDisplacement();
        final DenseMatrix64F precision = branchStatistics.getRawPrecision();

        if (!isFinite(actualization) || !isFinite(displacement) || !isFinite(precision)) {
            throw new IllegalStateException(
                    "Non-finite branch transition statistics in fillCanonicalTransition"
                            + "; actualizationFinite=" + isFinite(actualization)
                            + "; displacementFinite=" + isFinite(displacement)
                            + "; precisionFinite=" + isFinite(precision));
        }

        actualizationMatrix.set(actualization);
        displacementVector.set(displacement);
        branchPrecisionMatrix.set(precision);
        canonicalizeBranchPrecisionCovariancePair("fillCanonicalTransition");
        OUCanonicalTransitionBuilder.fillFromPrecisionMoments(
                actualizationMatrix, displacementVector, branchPrecisionMatrix, transition);
    }

    private void fillCanonicalTransitionFromKernel(final double branchLength,
                                                   final double[] optimum) {
        final double[][] transitionMatrix = transitionMatrixScratch;
        final double[][] transitionCovariance = transitionCovarianceScratch;
        final double[] transitionOffset = transitionOffsetScratch;
        fillTransitionMomentsFromKernel(branchLength, optimum, transitionMatrix, transitionOffset, transitionCovariance);

        // Keep branch moment buffers in sync with kernel-generated transition moments.
        // Several local-factor paths (e.g. recoverParentAbovePrecision) read these buffers.
        for (int i = 0; i < dimension; ++i) {
            displacementVector.unsafe_set(i, 0, transitionOffset[i]);
            for (int j = 0; j < dimension; ++j) {
                actualizationMatrix.unsafe_set(i, j, transitionMatrix[i][j]);
                branchCovarianceMatrix.unsafe_set(i, j, transitionCovariance[i][j]);
            }
        }
        canonicalizeBranchCovariancePrecisionPair(
                "fillCanonicalTransitionFromKernel(branchLength=" + branchLength + ")");

        if (branchLength <= 0.0) {
            // Zero-length branches can have singular covariance; use the local jittered path.
            OUCanonicalTransitionBuilder.fillFromMoments(
                    transitionMatrix,
                    transitionOffset,
                    transitionCovariance,
                    transition,
                    reducedCovarianceScratch,
                    reducedCholeskyScratch,
                    reducedPrecisionScratch,
                    reducedLowerInverseScratch,
                    reducedPrecisionScratch,
                    currentCovarianceScratch);
        } else {
            // Build canonical transition from the canonicalized precision pair to keep
            // all downstream consumers (local factors and parent-recovery algebra)
            // numerically consistent.
            OUCanonicalTransitionBuilder.fillFromPrecisionMoments(
                    actualizationMatrix, displacementVector, branchPrecisionMatrix, transition);
        }
    }

    private static void fillTransitionOffset(final double[][] transitionMatrix,
                                             final double[] optimum,
                                             final double[] out) {
        final int d = out.length;
        for (int i = 0; i < d; ++i) {
            double sum = optimum[i];
            for (int j = 0; j < d; ++j) {
                sum -= transitionMatrix[i][j] * optimum[j];
            }
            out[i] = sum;
        }
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
        final DenseMatrix64F aboveMean = parentAboveResolver.getAboveParentMeanVector();
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
        final DenseMatrix64F aboveMean = parentAboveResolver.getAboveParentMeanVector();
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

        CanonicalGaussianUtils.fillMomentsFromCanonical(currentPosterior, currentMeanScratch, currentCovarianceScratch);

        for (int i = 0; i < dimension; ++i) {
            final double xi = currentMeanScratch[i];
            final double yi = observedChild.unsafe_get(i, 0);
            contribution.dLogL_dInformationX[i] = xi;
            contribution.dLogL_dInformationY[i] = yi;
            for (int j = 0; j < dimension; ++j) {
                final double xj = currentMeanScratch[j];
                final double yj = observedChild.unsafe_get(j, 0);
                final double exx = currentCovarianceScratch[i][j] + xi * xj;
                final int ij = i * dimension + j;
                contribution.dLogL_dPrecisionXX[ij] = -0.5 * exx;
                contribution.dLogL_dPrecisionXY[ij] = -0.5 * (xi * yj);
                contribution.dLogL_dPrecisionYX[ij] = -0.5 * (yi * xj);
                contribution.dLogL_dPrecisionYY[ij] = -0.5 * (yi * yj);
            }
        }
        contribution.dLogL_dLogNormalizer = -1.0;
    }

    private void fillContributionForPartiallyObservedTip(final NormalSufficientStatistics above,
                                                         final NormalSufficientStatistics aboveParent,
                                                         final NormalSufficientStatistics below,
                                                         final int observedCount) {
        final DenseMatrix64F abovePrecision = recoverOrUseParentAbovePrecision(above, aboveParent);
        final DenseMatrix64F aboveMean = parentAboveResolver.getAboveParentMeanVector();
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

        clearContribution();
        for (int i = 0; i < dimension; ++i) {
            contribution.dLogL_dInformationX[i] = reducedMean[i];
            contribution.dLogL_dInformationY[i] = exactObservationMask[i]
                    ? observedChild.unsafe_get(i, 0)
                    : reducedMean[reducedIndexByTrait[i]];
        }

        for (int i = 0; i < dimension; ++i) {
            final double xi = reducedMean[i];
            final int reducedI = reducedIndexByTrait[i];
            for (int j = 0; j < dimension; ++j) {
                final double xj = reducedMean[j];
                final int reducedJ = reducedIndexByTrait[j];
                final double yi = contribution.dLogL_dInformationY[i];
                final double yj = contribution.dLogL_dInformationY[j];

                final double exx = reducedCovariance[i][j] + xi * xj;
                final double exy = exactObservationMask[j]
                        ? xi * yj
                        : reducedCovariance[i][reducedJ] + xi * yj;
                final double eyx = exactObservationMask[i]
                        ? yi * xj
                        : reducedCovariance[reducedI][j] + yi * xj;
                final double eyy = (exactObservationMask[i] || exactObservationMask[j])
                        ? yi * yj
                        : reducedCovariance[reducedI][reducedJ] + yi * yj;

                final int ij = i * dimension + j;
                contribution.dLogL_dPrecisionXX[ij] = -0.5 * exx;
                contribution.dLogL_dPrecisionXY[ij] = -0.5 * exy;
                contribution.dLogL_dPrecisionYX[ij] = -0.5 * eyx;
                contribution.dLogL_dPrecisionYY[ij] = -0.5 * eyy;
            }
        }
        contribution.dLogL_dLogNormalizer = -1.0;
    }

    private void fillContributionFromPairPosterior() {
        final int pairDimension = pairPosterior.getDimension();
        final int d = pairDimension / 2;
        fillMomentsFromCanonicalLocally(pairPosterior, reducedMeanScratch, reducedCovarianceScratch);

        for (int i = 0; i < d; ++i) {
            final double currentMeanI = reducedMeanScratch[i];
            final double nextMeanI = reducedMeanScratch[d + i];
            contribution.dLogL_dInformationX[i] = currentMeanI;
            contribution.dLogL_dInformationY[i] = nextMeanI;
            for (int j = 0; j < d; ++j) {
                final double currentMeanJ = reducedMeanScratch[j];
                final double nextMeanJ = reducedMeanScratch[d + j];
                final double currentSecondMoment =
                        reducedCovarianceScratch[i][j] + currentMeanI * currentMeanJ;
                final double crossSecondMoment =
                        reducedCovarianceScratch[d + i][j] + nextMeanI * currentMeanJ;
                final double nextSecondMoment =
                        reducedCovarianceScratch[d + i][d + j] + nextMeanI * nextMeanJ;
                final int ij = i * d + j;
                contribution.dLogL_dPrecisionXX[ij] = -0.5 * currentSecondMoment;
                contribution.dLogL_dPrecisionXY[ij] =
                        -0.5 * (reducedCovarianceScratch[i][d + j] + currentMeanI * nextMeanJ);
                contribution.dLogL_dPrecisionYX[ij] = -0.5 * crossSecondMoment;
                contribution.dLogL_dPrecisionYY[ij] = -0.5 * nextSecondMoment;
            }
        }
        contribution.dLogL_dLogNormalizer = -1.0;
    }

    private DenseMatrix64F recoverOrUseParentAbovePrecision(final NormalSufficientStatistics aboveChild,
                                                            final NormalSufficientStatistics aboveParent) {
        return parentAboveResolver.recoverOrUseParentAbovePrecision(aboveChild, aboveParent);
    }

    private DenseMatrix64F recoverParentAbovePrecision(final NormalSufficientStatistics above) {
        return parentAboveResolver.recoverParentAbovePrecision(above);
    }

    private void safeInvertSymmetricPositiveDefinite(final DenseMatrix64F source,
                                                     final DenseMatrix64F inverseOut) {
        safeInvertSymmetricPositiveDefinite(source, inverseOut, "unspecified");
    }

    private void canonicalizeBranchPrecisionCovariancePair(final String context) {
        symmetrizeInPlace(branchPrecisionMatrix);
        safeInvertSymmetricPositiveDefinite(
                branchPrecisionMatrix,
                branchCovarianceMatrix,
                context + ":branchPrecision->branchCovariance");
        safeInvertSymmetricPositiveDefinite(
                branchCovarianceMatrix,
                branchPrecisionMatrix,
                context + ":branchCovariance->branchPrecision");
        symmetrizeInPlace(branchPrecisionMatrix);
        symmetrizeInPlace(branchCovarianceMatrix);
    }

    private void canonicalizeBranchCovariancePrecisionPair(final String context) {
        symmetrizeInPlace(branchCovarianceMatrix);
        safeInvertSymmetricPositiveDefinite(
                branchCovarianceMatrix,
                branchPrecisionMatrix,
                context + ":branchCovariance->branchPrecision");
        safeInvertSymmetricPositiveDefinite(
                branchPrecisionMatrix,
                branchCovarianceMatrix,
                context + ":branchPrecision->branchCovariance");
        symmetrizeInPlace(branchPrecisionMatrix);
        symmetrizeInPlace(branchCovarianceMatrix);
    }

    private void safeInvertSymmetricPositiveDefinite(final DenseMatrix64F source,
                                                     final DenseMatrix64F inverseOut,
                                                     final String context) {
        CanonicalNumerics.safeInvertSymmetricPositiveDefinite(
                source,
                inverseOut,
                scratchMatrix,
                secondaryScratchMatrix,
                reducedPrecisionScratch,
                reducedCovarianceScratch,
                reducedCholeskyScratch,
                reducedLowerInverseScratch,
                NUMERICS_OPTIONS,
                context,
                parentAboveResolver.getSpdFailureDump());
    }

    private static void symmetrizeInPlace(final DenseMatrix64F matrix) {
        CanonicalNumerics.symmetrizeInPlace(matrix);
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
        final DenseMatrix64F precision = below.getRawPrecision();
        int observedCount = 0;
        for (int i = 0; i < dimension; ++i) {
            final double diagonal = precision.unsafe_get(i, i);
            if (Double.isInfinite(diagonal) && diagonal > 0.0) {
                exactObservationMask[i] = true;
                observedCount++;
            } else if (diagonal == 0.0) {
                exactObservationMask[i] = false;
            } else {
                return -1;
            }
            for (int j = 0; j < dimension; ++j) {
                if (i == j) {
                    continue;
                }
                if (precision.unsafe_get(i, j) != 0.0) {
                    return -1;
                }
            }
        }
        return observedCount;
    }

    private int collectObservationPartition(final int observedCount) {
        int observedCursor = 0;
        int missingCursor = 0;
        for (int i = 0; i < dimension; ++i) {
            if (exactObservationMask[i]) {
                observedIndexScratch[observedCursor++] = i;
                reducedIndexByTrait[i] = -1;
            } else {
                missingIndexScratch[missingCursor] = i;
                reducedIndexByTrait[i] = dimension + missingCursor;
                missingCursor++;
            }
        }
        if (observedCursor != observedCount) {
            throw new IllegalStateException("Observation partition count mismatch");
        }
        return missingCursor;
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
                final int observedIndex = observedIndexScratch[observed];
                information -= transition.precisionXY[iOffset + observedIndex] * observedChild.unsafe_get(observedIndex, 0);
            }
            reducedInformation[i] = information;
            for (int missing = 0; missing < missingCount; ++missing) {
                reducedPrecision[i][dimension + missing] =
                        transition.precisionXY[iOffset + missingIndexScratch[missing]];
            }
        }

        for (int missing = 0; missing < missingCount; ++missing) {
            final int childIndex = missingIndexScratch[missing];
            final int row = dimension + missing;
            double information = transition.informationY[childIndex];
            final int childOffset = childIndex * dimension;
            for (int observed = 0; observed < dimension - missingCount; ++observed) {
                final int observedIndex = observedIndexScratch[observed];
                information -= transition.precisionYY[childOffset + observedIndex] * observedChild.unsafe_get(observedIndex, 0);
            }
            reducedInformation[row] = information;
            for (int j = 0; j < dimension; ++j) {
                reducedPrecision[row][j] = transition.precisionYX[childOffset + j];
            }
            for (int otherMissing = 0; otherMissing < missingCount; ++otherMissing) {
                reducedPrecision[row][dimension + otherMissing] =
                        transition.precisionYY[childOffset + missingIndexScratch[otherMissing]];
            }
        }
    }

    private double evaluateFrozenPartiallyObservedLocalLogFactor(final DenseMatrix64F parentAbovePrecision,
                                                                 final DenseMatrix64F observedChild,
                                                                 final int observedCount) {
        final int missingCount = collectObservationPartition(observedCount);
        final int reducedDimension = dimension + missingCount;
        final double[][] reducedPrecision = reducedPrecisionScratch;
        final double[] reducedInformation = reducedInformationScratch;
        fillReducedCanonicalSystem(parentAbovePrecision, observedChild, missingCount, reducedPrecision, reducedInformation);

        double observedQuadratic = 0.0;
        double observedLinear = 0.0;
        for (int oi = 0; oi < observedCount; ++oi) {
            final int i = observedIndexScratch[oi];
            final double yi = observedChild.unsafe_get(i, 0);
            observedLinear += transition.informationY[i] * yi;
            final int iOffset = i * dimension;
            for (int oj = 0; oj < observedCount; ++oj) {
                final int j = observedIndexScratch[oj];
                observedQuadratic += yi * transition.precisionYY[iOffset + j] * observedChild.unsafe_get(j, 0);
            }
        }

        return normalizedLogNormalizer(reducedPrecision, reducedInformation, reducedDimension)
                - transition.logNormalizer
                - 0.5 * observedQuadratic
                + observedLinear;
    }

    private void clearContribution() {
        zero(contribution.dLogL_dInformationX);
        zero(contribution.dLogL_dInformationY);
        zero(contribution.dLogL_dPrecisionXX);
        zero(contribution.dLogL_dPrecisionXY);
        zero(contribution.dLogL_dPrecisionYX);
        zero(contribution.dLogL_dPrecisionYY);
        contribution.dLogL_dLogNormalizer = 0.0;
    }

    private double normalizedLogNormalizer(final double[][] precision,
                                           final double[] information,
                                           final int dimensionUsed) {
        final double[][] inverse = reducedCovarianceScratch;
        final double[] mean = reducedMeanScratch;
        final double logDet = invertSymmetricPositiveDefinite(precision, dimensionUsed, inverse);
        for (int i = 0; i < dimensionUsed; ++i) {
            double sum = 0.0;
            for (int j = 0; j < dimensionUsed; ++j) {
                sum += inverse[i][j] * information[j];
            }
            mean[i] = sum;
        }
        return 0.5 * (dimensionUsed * Math.log(2.0 * Math.PI) - logDet + dot(information, mean, dimensionUsed));
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

    private static double quadraticForm(final double[] matrix, final DenseMatrix64F vector) {
        double sum = 0.0;
        final int dimensionUsed = vector.numRows;
        for (int i = 0; i < dimensionUsed; ++i) {
            final int iOffset = i * dimensionUsed;
            for (int j = 0; j < dimensionUsed; ++j) {
                sum += vector.unsafe_get(i, 0) * matrix[iOffset + j] * vector.unsafe_get(j, 0);
            }
        }
        return sum;
    }

    private static double quadraticForm(final double[][] matrix, final DenseMatrix64F vector) {
        double sum = 0.0;
        for (int i = 0; i < vector.numRows; ++i) {
            for (int j = 0; j < vector.numRows; ++j) {
                sum += vector.unsafe_get(i, 0) * matrix[i][j] * vector.unsafe_get(j, 0);
            }
        }
        return sum;
    }

    private static double dot(final double[] left, final DenseMatrix64F right) {
        double sum = 0.0;
        for (int i = 0; i < left.length; ++i) {
            sum += left[i] * right.unsafe_get(i, 0);
        }
        return sum;
    }

    private static double dot(final double[] left, final double[] right) {
        double sum = 0.0;
        for (int i = 0; i < left.length; ++i) {
            sum += left[i] * right[i];
        }
        return sum;
    }

    private static double dot(final double[] left, final double[] right, final int dimensionUsed) {
        double sum = 0.0;
        for (int i = 0; i < dimensionUsed; ++i) {
            sum += left[i] * right[i];
        }
        return sum;
    }

}

package dr.evomodel.treedatalikelihood.continuous;

import dr.evomodel.treedatalikelihood.preorder.BranchSufficientStatistics;
import dr.evomodel.treedatalikelihood.preorder.MatrixSufficientStatistics;
import dr.evomodel.treedatalikelihood.preorder.NormalSufficientStatistics;
import dr.inference.timeseries.engine.gaussian.CanonicalBranchMessageContribution;
import dr.inference.timeseries.engine.gaussian.CanonicalLocalTransitionAdjoints;
import dr.inference.timeseries.engine.gaussian.CanonicalTransitionAdjointUtils;
import dr.inference.timeseries.gaussian.OUProcessModel;
import dr.inference.timeseries.representation.CanonicalGaussianState;
import dr.inference.timeseries.representation.CanonicalGaussianTransition;
import dr.inference.timeseries.representation.CanonicalGaussianUtils;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;

/**
 * Self-contained local tree-to-time-series canonical wiring for one OU branch.
 *
 * <p>This class is intentionally narrow: given tree-side branch sufficient statistics,
 * it reconstructs the canonical branch factor that the time-series machinery expects,
 * handles the exact observed-tip conditioning limit, and produces local canonical
 * transition adjoints. That keeps the branch factor construction explicit and testable
 * outside the larger delegate/gradient classes.</p>
 */
public final class TimeSeriesOUCanonicalBranchWiring {
    private static final String DISABLE_EXACT_TIP_SHORTCUT_PROPERTY =
            "beast.experimental.disableExactTipShortcut";
    private static final String SPD_DEBUG_DUMP_PROPERTY =
            "beast.debug.ou.spdDump";
    private static final String NONFINITE_BRANCH_STATS_DEBUG_PROPERTY =
            "beast.debug.ou.nonfiniteBranchStats";
    private static final String FORCE_STRICT_SPD_INVERSION_PROPERTY =
            "beast.debug.ou.forceStrictSpdInverse";
    private static final String COMPARE_PARENT_ABOVE_RECOVERED_PROPERTY =
            "beast.debug.ou.compareParentAboveRecovered";
    private static final String COMPARE_PARENT_ABOVE_RECOVERED_FAIL_PROPERTY =
            "beast.debug.ou.compareParentAboveRecovered.failOnMismatch";
    private static final String COMPARE_PARENT_ABOVE_RECOVERED_TOLERANCE_PROPERTY =
            "beast.debug.ou.compareParentAboveRecovered.tolerance";
    private static final String COMPARE_PARENT_ABOVE_RECOVERED_MAX_REPORTS_PROPERTY =
            "beast.debug.ou.compareParentAboveRecovered.maxReports";
    private static final String FORCE_RECOVER_PARENT_ABOVE_FROM_CHILD_PROPERTY =
            "beast.debug.ou.forceRecoverParentAboveFromChild";
    private static final String ALLOW_INFINITE_PARENT_ABOVE_RECOVERY_PROPERTY =
            "beast.debug.ou.allowInfiniteParentAboveRecovery";

    private static final double SYMMETRIC_JITTER_RELATIVE = 1.0e-12;
    private static final double SYMMETRIC_JITTER_ABSOLUTE = 1.0e-12;

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
    private final DenseMatrix64F aboveChildPrecisionMatrix;
    private final DenseMatrix64F aboveChildCovarianceMatrix;
    private final DenseMatrix64F aboveParentPrecisionMatrix;
    private final DenseMatrix64F aboveParentCovarianceMatrix;
    private final DenseMatrix64F actualizationInverseMatrix;
    private final DenseMatrix64F aboveParentMeanVector;
    private final DenseMatrix64F centeredChildMeanVector;
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
    private int parentAboveCompareReportCount = 0;

    public TimeSeriesOUCanonicalBranchWiring(final TimeSeriesOUGaussianBranchTransitionProvider branchTransitionProvider) {
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
        this.aboveChildPrecisionMatrix = new DenseMatrix64F(dimension, dimension);
        this.aboveChildCovarianceMatrix = new DenseMatrix64F(dimension, dimension);
        this.aboveParentPrecisionMatrix = new DenseMatrix64F(dimension, dimension);
        this.aboveParentCovarianceMatrix = new DenseMatrix64F(dimension, dimension);
        this.actualizationInverseMatrix = new DenseMatrix64F(dimension, dimension);
        this.aboveParentMeanVector = new DenseMatrix64F(dimension, 1);
        this.centeredChildMeanVector = new DenseMatrix64F(dimension, 1);
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
        final DenseMatrix64F parentAboveMean = aboveParentMeanVector;
        return evaluateFrozenLocalLogFactorWithResolvedParentPrecision(parentAbovePrecision, parentAboveMean, below);
    }

    public double evaluateFrozenLocalLogFactor(final double branchLength,
                                               final double[] optimum,
                                               final NormalSufficientStatistics aboveChild,
                                               final NormalSufficientStatistics aboveParent,
                                               final NormalSufficientStatistics below) {
        fillCanonicalTransitionFromKernel(branchLength, optimum);
        final DenseMatrix64F parentAbovePrecision = recoverOrUseParentAbovePrecision(aboveChild, aboveParent);
        final DenseMatrix64F parentAboveMean = aboveParentMeanVector;
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

        for (int i = 0; i < dimension; ++i) {
            double infoY = 0.0;
            final int iOffset = i * dimension;
            for (int j = 0; j < dimension; ++j) {
                final double p = branchPrecisionMatrix.unsafe_get(i, j);
                transition.precisionYY[iOffset + j] = p;
                transition.precisionYX[iOffset + j] = -multiplyEntry(branchPrecisionMatrix, actualization, i, j);
                transition.precisionXY[j * dimension + i] = transition.precisionYX[iOffset + j];
                transition.precisionXX[iOffset + j] =
                        multiplyEntryTranspose(actualization, branchPrecisionMatrix, actualization, i, j);
                infoY += p * displacement.unsafe_get(j, 0);
            }
            transition.informationY[i] = infoY;
        }
        for (int i = 0; i < dimension; ++i) {
            double infoX = 0.0;
            for (int j = 0; j < dimension; ++j) {
                infoX -= actualization.unsafe_get(j, i) * transition.informationY[j];
            }
            transition.informationX[i] = infoX;
        }
        transition.logNormalizer = 0.0;
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
            fillTransitionFromMomentsLocally(
                    transitionMatrix, transitionOffset, transitionCovariance, transition);
        } else {
            // Build canonical transition from the canonicalized precision pair to keep
            // all downstream consumers (local factors and parent-recovery algebra)
            // numerically consistent.
            for (int i = 0; i < dimension; ++i) {
                double infoY = 0.0;
                final int iOffset = i * dimension;
                for (int j = 0; j < dimension; ++j) {
                    final double p = branchPrecisionMatrix.unsafe_get(i, j);
                    transition.precisionYY[iOffset + j] = p;
                    transition.precisionYX[iOffset + j] = -multiplyEntry(branchPrecisionMatrix, actualizationMatrix, i, j);
                    transition.precisionXY[j * dimension + i] = transition.precisionYX[iOffset + j];
                    transition.precisionXX[iOffset + j] =
                            multiplyEntryTranspose(actualizationMatrix, branchPrecisionMatrix, actualizationMatrix, i, j);
                    infoY += p * displacementVector.unsafe_get(j, 0);
                }
                transition.informationY[i] = infoY;
            }
            for (int i = 0; i < dimension; ++i) {
                double infoX = 0.0;
                for (int j = 0; j < dimension; ++j) {
                    infoX -= actualizationMatrix.unsafe_get(j, i) * transition.informationY[j];
                }
                transition.informationX[i] = infoX;
            }
            transition.logNormalizer = 0.0;
        }
    }

    private void fillTransitionFromMomentsLocally(final double[][] transitionMatrix,
                                                  final double[] transitionOffset,
                                                  final double[][] transitionCovariance,
                                                  final CanonicalGaussianTransition out) {
        final double[][] precision = reducedCovarianceScratch;
        final double[][] transitionTranspose = reducedPrecisionScratch;
        final double[][] tmp = currentCovarianceScratch;

        final double logDet = invertSymmetricPositiveDefinite(transitionCovariance, dimension, precision);
        transpose(transitionMatrix, transitionTranspose, dimension);

        multiply(transitionTranspose, precision, tmp, dimension);
        multiply(tmp, transitionMatrix, out.precisionXX, dimension);

        multiply(transitionTranspose, precision, out.precisionXY, dimension);
        scaleInPlace(out.precisionXY, -1.0, dimension * dimension);

        multiply(precision, transitionMatrix, out.precisionYX, dimension);
        scaleInPlace(out.precisionYX, -1.0, dimension * dimension);

        copyMatrix(precision, out.precisionYY, dimension);

        multiply(precision, transitionOffset, out.informationY, dimension);
        multiply(transitionTranspose, out.informationY, out.informationX, dimension);
        scaleInPlace(out.informationX, -1.0, dimension);

        out.logNormalizer =
                0.5 * (dimension * Math.log(2.0 * Math.PI)
                        + logDet
                        + dot(transitionOffset, out.informationY, dimension));
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
        for (double value : values) {
            if (!Double.isFinite(value)) {
                return false;
            }
        }
        return true;
    }

    private static boolean isFinite(final double[][] values) {
        for (double[] row : values) {
            if (!isFinite(row)) {
                return false;
            }
        }
        return true;
    }

    private static double maxAbs(final double[] values) {
        double max = 0.0;
        for (double value : values) {
            max = Math.max(max, Math.abs(value));
        }
        return max;
    }

    private static void transpose(final double[][] source,
                                  final double[][] target,
                                  final int dimensionUsed) {
        for (int i = 0; i < dimensionUsed; ++i) {
            for (int j = 0; j < dimensionUsed; ++j) {
                target[j][i] = source[i][j];
            }
        }
    }

    private static void multiply(final double[][] left,
                                 final double[][] right,
                                 final double[][] out,
                                 final int dimensionUsed) {
        for (int i = 0; i < dimensionUsed; ++i) {
            for (int j = 0; j < dimensionUsed; ++j) {
                double sum = 0.0;
                for (int k = 0; k < dimensionUsed; ++k) {
                    sum += left[i][k] * right[k][j];
                }
                out[i][j] = sum;
            }
        }
    }

    private static void multiply(final double[][] matrix,
                                 final double[] vector,
                                 final double[] out,
                                 final int dimensionUsed) {
        for (int i = 0; i < dimensionUsed; ++i) {
            double sum = 0.0;
            for (int j = 0; j < dimensionUsed; ++j) {
                sum += matrix[i][j] * vector[j];
            }
            out[i] = sum;
        }
    }

    private static void multiply(final double[][] left,
                                 final double[][] right,
                                 final double[] out,
                                 final int dimensionUsed) {
        for (int i = 0; i < dimensionUsed; ++i) {
            final int iOffset = i * dimensionUsed;
            for (int j = 0; j < dimensionUsed; ++j) {
                double sum = 0.0;
                for (int k = 0; k < dimensionUsed; ++k) {
                    sum += left[i][k] * right[k][j];
                }
                out[iOffset + j] = sum;
            }
        }
    }

    private static void copyMatrix(final double[][] source,
                                   final double[][] target,
                                   final int dimensionUsed) {
        for (int i = 0; i < dimensionUsed; ++i) {
            System.arraycopy(source[i], 0, target[i], 0, dimensionUsed);
        }
    }

    private static void copyMatrix(final double[][] source,
                                   final double[] target,
                                   final int dimensionUsed) {
        for (int i = 0; i < dimensionUsed; ++i) {
            System.arraycopy(source[i], 0, target, i * dimensionUsed, dimensionUsed);
        }
    }

    private static void scaleInPlace(final double[][] matrix,
                                     final double factor,
                                     final int dimensionUsed) {
        for (int i = 0; i < dimensionUsed; ++i) {
            for (int j = 0; j < dimensionUsed; ++j) {
                matrix[i][j] *= factor;
            }
        }
    }

    private static void scaleInPlace(final double[] vector,
                                     final double factor,
                                     final int dimensionUsed) {
        for (int i = 0; i < dimensionUsed; ++i) {
            vector[i] *= factor;
        }
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
        final DenseMatrix64F aboveMean = aboveParentMeanVector;
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
        final DenseMatrix64F aboveMean = aboveParentMeanVector;
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
        final DenseMatrix64F aboveMean = aboveParentMeanVector;
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
        if (Boolean.getBoolean(FORCE_RECOVER_PARENT_ABOVE_FROM_CHILD_PROPERTY)) {
            return recoverParentAbovePrecision(aboveChild);
        }
        if (aboveParent == null) {
            throw new IllegalStateException(
                    "Missing parent-above message for canonical branch contribution."
                            + " childAbovePrecisionFinite=" + isFinite(aboveChild.getRawPrecision())
                            + " childAboveMeanFinite=" + isFinite(aboveChild.getRawMean()));
        }
        if (Boolean.getBoolean(COMPARE_PARENT_ABOVE_RECOVERED_PROPERTY)) {
            compareParentAboveAgainstRecovered(aboveChild, aboveParent);
        }
        final DenseMatrix64F parentPrecisionRaw = aboveParent.getRawPrecision();
        final DenseMatrix64F parentMeanRaw = aboveParent.getRawMean();
        final boolean parentMeanFinite = isFinite(parentMeanRaw);
        final boolean parentPrecisionHasNaN = hasNaN(parentPrecisionRaw);
        final boolean parentPrecisionHasInfinity = hasInfinity(parentPrecisionRaw);
        if (!parentMeanFinite || parentPrecisionHasNaN) {
            throw new IllegalStateException(
                    "Invalid parent-above message."
                            + " parentMeanFinite=" + parentMeanFinite
                            + " parentPrecisionHasNaN=" + parentPrecisionHasNaN
                            + " parentPrecisionHasInfinity=" + parentPrecisionHasInfinity
                            + " parentPrecisionSummary=" + summarizeDenseMatrix(parentPrecisionRaw)
                            + " parentMeanSummary=" + summarizeDenseMatrix(parentMeanRaw)
                            + " childAbovePrecisionFinite=" + isFinite(aboveChild.getRawPrecision())
                            + " childAboveMeanFinite=" + isFinite(aboveChild.getRawMean())
                            + " childAbovePrecisionSummary=" + summarizeDenseMatrix(aboveChild.getRawPrecision())
                            + " childAboveMeanSummary=" + summarizeDenseMatrix(aboveChild.getRawMean()));
        }
        if (parentPrecisionHasInfinity) {
            if (!Boolean.parseBoolean(System.getProperty(ALLOW_INFINITE_PARENT_ABOVE_RECOVERY_PROPERTY, "true"))) {
                throw new IllegalStateException(
                        "Parent-above precision contains infinite entries and infinite-recovery is disabled."
                                + " parentPrecisionSummary=" + summarizeDenseMatrix(parentPrecisionRaw)
                                + " parentMeanSummary=" + summarizeDenseMatrix(parentMeanRaw));
            }
            // Parent-above messages can be represented as exact constraints (infinite precision)
            // in mixed missing-data paths. Convert to a finite equivalent recovered from child-above
            // for stable downstream canonical algebra. Keep the supplied parent mean fixed:
            // recovering mean from child-above introduces parameter dependence that is not
            // part of the original parent-above message.
            final DenseMatrix64F recoveredPrecision = recoverParentAbovePrecision(aboveChild);
            aboveParentMeanVector.set(parentMeanRaw);
            return recoveredPrecision;
        }
        aboveParentPrecisionMatrix.set(parentPrecisionRaw);
        aboveParentMeanVector.set(parentMeanRaw);
        return aboveParentPrecisionMatrix;
    }

    private void compareParentAboveAgainstRecovered(final NormalSufficientStatistics aboveChild,
                                                    final NormalSufficientStatistics aboveParent) {
        final DenseMatrix64F suppliedPrecision = aboveParent.getRawPrecision();
        final DenseMatrix64F suppliedMean = aboveParent.getRawMean();

        final DenseMatrix64F recoveredPrecision = recoverParentAbovePrecision(aboveChild);
        final DenseMatrix64F recoveredMean = aboveParentMeanVector;

        final MatrixComparison precisionComparison = compareMatricesAllowingInf(suppliedPrecision, recoveredPrecision);
        final MatrixComparison meanComparison = compareMatricesAllowingInf(suppliedMean, recoveredMean);

        final double tolerance = Double.parseDouble(
                System.getProperty(COMPARE_PARENT_ABOVE_RECOVERED_TOLERANCE_PROPERTY, "1e-6"));
        final boolean mismatch =
                precisionComparison.infPatternMismatchCount > 0
                        || precisionComparison.nanPatternMismatchCount > 0
                        || precisionComparison.maxFiniteAbsDiff > tolerance
                        || meanComparison.infPatternMismatchCount > 0
                        || meanComparison.nanPatternMismatchCount > 0
                        || meanComparison.maxFiniteAbsDiff > tolerance;

        final int maxReports = Integer.parseInt(
                System.getProperty(COMPARE_PARENT_ABOVE_RECOVERED_MAX_REPORTS_PROPERTY, "20"));
        if (mismatch || parentAboveCompareReportCount < maxReports) {
            if (parentAboveCompareReportCount < maxReports) {
                parentAboveCompareReportCount++;
                System.err.println(
                        "parentAboveCompare"
                                + " call=" + parentAboveCompareReportCount
                                + " precision=" + precisionComparison
                                + " mean=" + meanComparison
                                + " tolerance=" + tolerance);
            }
        }

        if (mismatch && Boolean.getBoolean(COMPARE_PARENT_ABOVE_RECOVERED_FAIL_PROPERTY)) {
            throw new IllegalStateException(
                    "parent-above supplied vs recovered mismatch"
                            + " precision=" + precisionComparison
                            + " mean=" + meanComparison
                            + " tolerance=" + tolerance);
        }
    }

    private static MatrixComparison compareMatricesAllowingInf(final DenseMatrix64F supplied,
                                                               final DenseMatrix64F recovered) {
        final double[] suppliedData = supplied.getData();
        final double[] recoveredData = recovered.getData();
        final int length = suppliedData.length;

        int finitePairCount = 0;
        int infPatternMismatchCount = 0;
        int nanPatternMismatchCount = 0;
        double maxFiniteAbsDiff = 0.0;

        for (int i = 0; i < length; ++i) {
            final double a = suppliedData[i];
            final double b = recoveredData[i];

            final boolean aNaN = Double.isNaN(a);
            final boolean bNaN = Double.isNaN(b);
            if (aNaN || bNaN) {
                if (!(aNaN && bNaN)) {
                    nanPatternMismatchCount++;
                }
                continue;
            }

            final boolean aInf = Double.isInfinite(a);
            final boolean bInf = Double.isInfinite(b);
            if (aInf || bInf) {
                if (!(aInf && bInf && Math.signum(a) == Math.signum(b))) {
                    infPatternMismatchCount++;
                }
                continue;
            }

            finitePairCount++;
            maxFiniteAbsDiff = Math.max(maxFiniteAbsDiff, Math.abs(a - b));
        }

        return new MatrixComparison(
                maxFiniteAbsDiff,
                finitePairCount,
                infPatternMismatchCount,
                nanPatternMismatchCount);
    }

    private static final class MatrixComparison {
        private final double maxFiniteAbsDiff;
        private final int finitePairCount;
        private final int infPatternMismatchCount;
        private final int nanPatternMismatchCount;

        private MatrixComparison(final double maxFiniteAbsDiff,
                                 final int finitePairCount,
                                 final int infPatternMismatchCount,
                                 final int nanPatternMismatchCount) {
            this.maxFiniteAbsDiff = maxFiniteAbsDiff;
            this.finitePairCount = finitePairCount;
            this.infPatternMismatchCount = infPatternMismatchCount;
            this.nanPatternMismatchCount = nanPatternMismatchCount;
        }

        @Override
        public String toString() {
            return "{maxFiniteAbsDiff=" + maxFiniteAbsDiff
                    + ",finitePairCount=" + finitePairCount
                    + ",infPatternMismatchCount=" + infPatternMismatchCount
                    + ",nanPatternMismatchCount=" + nanPatternMismatchCount
                    + "}";
        }
    }

    private DenseMatrix64F recoverParentAbovePrecision(final NormalSufficientStatistics above) {
        aboveChildPrecisionMatrix.set(above.getRawPrecision());
        safeInvertSymmetricPositiveDefinite(
                aboveChildPrecisionMatrix,
                aboveChildCovarianceMatrix,
                "recoverParentAbovePrecision:aboveChildPrecision->aboveChildCovariance");

        for (int i = 0; i < dimension; ++i) {
            centeredChildMeanVector.unsafe_set(i, 0,
                    above.getRawMean().unsafe_get(i, 0) - displacementVector.unsafe_get(i, 0));
        }

        // Keep this path smooth and deterministic:
        // always regularize/invert the same stabilized actualization matrix
        // instead of branching between solve/fallback paths.
        secondaryScratchMatrix.set(actualizationMatrix);
        final double solveJitter = Math.max(
                SYMMETRIC_JITTER_ABSOLUTE,
                SYMMETRIC_JITTER_RELATIVE * Math.max(1.0, maxAbsDiagonal(secondaryScratchMatrix)));
        addDiagonalJitter(secondaryScratchMatrix, solveJitter);
        safeInvert(secondaryScratchMatrix, actualizationInverseMatrix);
        CommonOps.mult(actualizationInverseMatrix, centeredChildMeanVector, aboveParentMeanVector);

        aboveParentCovarianceMatrix.set(aboveChildCovarianceMatrix);
        CommonOps.subtractEquals(aboveParentCovarianceMatrix, branchCovarianceMatrix);
        symmetrizeInPlace(aboveParentCovarianceMatrix);
        CommonOps.mult(actualizationInverseMatrix, aboveParentCovarianceMatrix, scratchMatrix);
        CommonOps.multTransB(scratchMatrix, actualizationInverseMatrix, aboveParentCovarianceMatrix);
        symmetrizeInPlace(aboveParentCovarianceMatrix);
        safeInvertSymmetricPositiveDefinite(
                aboveParentCovarianceMatrix,
                aboveParentPrecisionMatrix,
                "recoverParentAbovePrecision:aboveParentCovariance->aboveParentPrecision");
        return aboveParentPrecisionMatrix;
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
        symmetrizeCopy(source, scratchMatrix);
        final double jitterBase = Math.max(
                SYMMETRIC_JITTER_ABSOLUTE,
                SYMMETRIC_JITTER_RELATIVE * Math.max(1.0, maxAbsDiagonal(scratchMatrix)));

        if (Boolean.getBoolean(FORCE_STRICT_SPD_INVERSION_PROPERTY)) {
            if (invertSymmetricPositiveDefiniteStrictFallback(scratchMatrix, inverseOut, jitterBase)) {
                return;
            }
            emitSpdFailureDebugDump(context, source, scratchMatrix, jitterBase);
            throw new IllegalStateException(
                    "Failed strict SPD inversion with fallback; context=" + context);
        }

        double jitter = 0.0;
        for (int attempt = 0; attempt < 8; ++attempt) {
            secondaryScratchMatrix.set(scratchMatrix);
            if (jitter > 0.0) {
                addDiagonalJitter(secondaryScratchMatrix, jitter);
            }
            if (CommonOps.invert(secondaryScratchMatrix, inverseOut) && isFinite(inverseOut)) {
                symmetrizeInPlace(inverseOut);
                return;
            }
            jitter = jitter == 0.0 ? jitterBase : 10.0 * jitter;
        }
        if (invertSymmetricPositiveDefiniteStrictFallback(scratchMatrix, inverseOut, jitterBase)) {
            return;
        }
        emitSpdFailureDebugDump(context, source, scratchMatrix, jitterBase);
        throw new IllegalStateException(
                "Failed to invert symmetric positive definite matrix stably; context=" + context);
    }

    private boolean invertSymmetricPositiveDefiniteStrictFallback(final DenseMatrix64F symmetricSource,
                                                                  final DenseMatrix64F inverseOut,
                                                                  final double jitterBase) {
        final int d = symmetricSource.numRows;
        final double[][] symmetric = reducedPrecisionScratch;
        final double[][] adjusted = reducedCovarianceScratch;
        final double[][] cholesky = reducedCholeskyScratch;
        final double[][] lowerInverse = reducedLowerInverseScratch;

        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d; ++j) {
                symmetric[i][j] = symmetricSource.unsafe_get(i, j);
            }
        }

        double jitter = 0.0;
        for (int attempt = 0; attempt < 12; ++attempt) {
            copySquare(symmetric, adjusted, d);
            if (jitter > 0.0) {
                for (int i = 0; i < d; ++i) {
                    adjusted[i][i] += jitter;
                }
            }
            if (invertSymmetricPositiveDefiniteStrict(adjusted, cholesky, lowerInverse, inverseOut, d)) {
                symmetrizeInPlace(inverseOut);
                return true;
            }
            jitter = jitter == 0.0 ? jitterBase : 10.0 * jitter;
        }
        return false;
    }

    private static boolean invertSymmetricPositiveDefiniteStrict(final double[][] matrix,
                                                                 final double[][] cholesky,
                                                                 final double[][] lowerInverse,
                                                                 final DenseMatrix64F inverseOut,
                                                                 final int dimensionUsed) {
        final int d = dimensionUsed;
        copySquare(matrix, cholesky, d);

        for (int i = 0; i < d; ++i) {
            for (int j = 0; j <= i; ++j) {
                double sum = cholesky[i][j];
                for (int k = 0; k < j; ++k) {
                    sum -= cholesky[i][k] * cholesky[j][k];
                }
                if (i == j) {
                    if (!(sum > 0.0) || !Double.isFinite(sum)) {
                        return false;
                    }
                    cholesky[i][j] = Math.sqrt(sum);
                } else {
                    cholesky[i][j] = sum / cholesky[j][j];
                }
            }
            for (int j = i + 1; j < d; ++j) {
                cholesky[i][j] = 0.0;
            }
        }

        for (int column = 0; column < d; ++column) {
            for (int row = 0; row < d; ++row) {
                double sum = row == column ? 1.0 : 0.0;
                for (int k = 0; k < row; ++k) {
                    sum -= cholesky[row][k] * lowerInverse[k][column];
                }
                lowerInverse[row][column] = sum / cholesky[row][row];
            }
        }

        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d; ++j) {
                double sum = 0.0;
                for (int k = 0; k < d; ++k) {
                    sum += lowerInverse[k][i] * lowerInverse[k][j];
                }
                inverseOut.unsafe_set(i, j, sum);
            }
        }
        return isFinite(inverseOut);
    }

    private void safeInvert(final DenseMatrix64F source,
                            final DenseMatrix64F inverseOut) {
        final double jitterBase = Math.max(
                SYMMETRIC_JITTER_ABSOLUTE,
                SYMMETRIC_JITTER_RELATIVE * Math.max(1.0, maxAbsDiagonal(source)));
        double jitter = jitterBase;
        for (int attempt = 0; attempt < 8; ++attempt) {
            secondaryScratchMatrix.set(source);
            addDiagonalJitter(secondaryScratchMatrix, jitter);
            if (CommonOps.invert(secondaryScratchMatrix, inverseOut) && isFinite(inverseOut)) {
                return;
            }
            jitter *= 10.0;
        }
        throw new IllegalStateException("Failed to invert matrix stably");
    }

    private static void symmetrizeCopy(final DenseMatrix64F source,
                                       final DenseMatrix64F destination) {
        final int d = source.numRows;
        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d; ++j) {
                destination.unsafe_set(i, j,
                        0.5 * (source.unsafe_get(i, j) + source.unsafe_get(j, i)));
            }
        }
    }

    private static void symmetrizeInPlace(final DenseMatrix64F matrix) {
        final int d = matrix.numRows;
        for (int i = 0; i < d; ++i) {
            for (int j = i + 1; j < d; ++j) {
                final double value = 0.5 * (matrix.unsafe_get(i, j) + matrix.unsafe_get(j, i));
                matrix.unsafe_set(i, j, value);
                matrix.unsafe_set(j, i, value);
            }
        }
    }

    private static void addDiagonalJitter(final DenseMatrix64F matrix,
                                          final double jitter) {
        for (int i = 0; i < matrix.numRows; ++i) {
            matrix.unsafe_set(i, i, matrix.unsafe_get(i, i) + jitter);
        }
    }

    private static double maxAbsDiagonal(final DenseMatrix64F matrix) {
        double max = 0.0;
        for (int i = 0; i < matrix.numRows; ++i) {
            max = Math.max(max, Math.abs(matrix.unsafe_get(i, i)));
        }
        return max;
    }

    private static boolean isFinite(final DenseMatrix64F matrix) {
        final double[] data = matrix.getData();
        for (double value : data) {
            if (!Double.isFinite(value)) {
                return false;
            }
        }
        return true;
    }

    private static boolean hasNaN(final DenseMatrix64F matrix) {
        final double[] data = matrix.getData();
        for (double value : data) {
            if (Double.isNaN(value)) {
                return true;
            }
        }
        return false;
    }

    private static boolean hasInfinity(final DenseMatrix64F matrix) {
        final double[] data = matrix.getData();
        for (double value : data) {
            if (Double.isInfinite(value)) {
                return true;
            }
        }
        return false;
    }

    private static String summarizeDenseMatrix(final DenseMatrix64F matrix) {
        int nanCount = 0;
        int posInfCount = 0;
        int negInfCount = 0;
        double minFinite = Double.POSITIVE_INFINITY;
        double maxFinite = Double.NEGATIVE_INFINITY;
        final double[] data = matrix.getData();
        for (double value : data) {
            if (Double.isNaN(value)) {
                nanCount++;
            } else if (value == Double.POSITIVE_INFINITY) {
                posInfCount++;
            } else if (value == Double.NEGATIVE_INFINITY) {
                negInfCount++;
            } else {
                minFinite = Math.min(minFinite, value);
                maxFinite = Math.max(maxFinite, value);
            }
        }
        if (minFinite == Double.POSITIVE_INFINITY) {
            minFinite = Double.NaN;
            maxFinite = Double.NaN;
        }
        return "{rows=" + matrix.numRows
                + ",cols=" + matrix.numCols
                + ",nan=" + nanCount
                + ",posInf=" + posInfCount
                + ",negInf=" + negInfCount
                + ",minFinite=" + minFinite
                + ",maxFinite=" + maxFinite
                + "}";
    }

    private static void copySquare(final double[][] source,
                                   final double[][] target,
                                   final int dimensionUsed) {
        for (int i = 0; i < dimensionUsed; ++i) {
            System.arraycopy(source[i], 0, target[i], 0, dimensionUsed);
        }
    }

    private static void symmetrizeInPlace(final double[][] matrix,
                                          final int dimensionUsed) {
        for (int i = 0; i < dimensionUsed; ++i) {
            for (int j = i + 1; j < dimensionUsed; ++j) {
                final double value = 0.5 * (matrix[i][j] + matrix[j][i]);
                matrix[i][j] = value;
                matrix[j][i] = value;
            }
        }
    }

    private static boolean isFinite(final double[][] matrix,
                                    final int dimensionUsed) {
        for (int i = 0; i < dimensionUsed; ++i) {
            for (int j = 0; j < dimensionUsed; ++j) {
                if (!Double.isFinite(matrix[i][j])) {
                    return false;
                }
            }
        }
        return true;
    }

    private static double maxAbsDiagonal(final double[][] matrix,
                                         final int dimensionUsed) {
        double max = 0.0;
        for (int i = 0; i < dimensionUsed; ++i) {
            max = Math.max(max, Math.abs(matrix[i][i]));
        }
        return max;
    }

    private static double minDiagonal(final double[][] matrix,
                                      final int dimensionUsed) {
        double min = Double.POSITIVE_INFINITY;
        for (int i = 0; i < dimensionUsed; ++i) {
            min = Math.min(min, matrix[i][i]);
        }
        return min;
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

    private static double multiplyEntry(final DenseMatrix64F left,
                                        final DenseMatrix64F right,
                                        final int row,
                                        final int column) {
        double sum = 0.0;
        final int inner = left.numCols;
        for (int k = 0; k < inner; ++k) {
            sum += left.unsafe_get(row, k) * right.unsafe_get(k, column);
        }
        return sum;
    }

    private static double multiplyEntryTranspose(final DenseMatrix64F left,
                                                 final DenseMatrix64F middle,
                                                 final DenseMatrix64F right,
                                                 final int row,
                                                 final int column) {
        double sum = 0.0;
        final int inner = middle.numRows;
        for (int i = 0; i < inner; ++i) {
            double leftMiddle = 0.0;
            for (int k = 0; k < inner; ++k) {
                leftMiddle += left.unsafe_get(k, row) * middle.unsafe_get(k, i);
            }
            sum += leftMiddle * right.unsafe_get(i, column);
        }
        return sum;
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
        // Keep source symmetrization and inversion output on disjoint buffers.
        // Some call sites pass `reducedCovarianceScratch` as inverseOut; using it as
        // the symmetrized source can corrupt retry attempts after a failed factorization.
        final double[][] symmetric = reducedCholeskyScratch;
        final double[][] cholesky = reducedPrecisionScratch;
        final double[][] lowerInverse = reducedLowerInverseScratch;

        for (int i = 0; i < dimensionUsed; ++i) {
            for (int j = 0; j < dimensionUsed; ++j) {
                symmetric[i][j] = 0.5 * (matrix[i][j] + matrix[j][i]);
            }
        }

        final double jitterBase = Math.max(
                SYMMETRIC_JITTER_ABSOLUTE,
                SYMMETRIC_JITTER_RELATIVE * Math.max(1.0, maxAbsDiagonal(symmetric, dimensionUsed)));

        final double lowerBound = gershgorinLowerBound(symmetric, dimensionUsed);
        double jitter = 0.0;
        for (int attempt = 0; attempt < 12; ++attempt) {
            for (int i = 0; i < dimensionUsed; ++i) {
                System.arraycopy(symmetric[i], 0, cholesky[i], 0, dimensionUsed);
            }
            if (jitter > 0.0) {
                for (int i = 0; i < dimensionUsed; ++i) {
                    cholesky[i][i] += jitter;
                }
            }
            final double logDet = invertSymmetricPositiveDefiniteStrict(cholesky, lowerInverse, inverseOut, dimensionUsed);
            if (Double.isFinite(logDet) && isFinite(inverseOut, dimensionUsed)) {
                symmetrizeInPlace(inverseOut, dimensionUsed);
                return logDet;
            }
            if (jitter == 0.0) {
                jitter = lowerBound > 0.0 ? jitterBase : (-lowerBound + jitterBase);
            } else {
                jitter *= 10.0;
            }
        }
        throw new IllegalStateException(
                "Failed to invert symmetric positive definite matrix stably"
                        + "; dim=" + dimensionUsed
                        + "; finite=" + isFinite(symmetric, dimensionUsed)
                        + "; gershgorinLowerBound=" + gershgorinLowerBound(symmetric, dimensionUsed)
                        + "; minDiag=" + minDiagonal(symmetric, dimensionUsed)
                        + "; maxAbsDiag=" + maxAbsDiagonal(symmetric, dimensionUsed));
    }

    private double invertSymmetricPositiveDefiniteStrict(final double[][] cholesky,
                                                         final double[][] lowerInverse,
                                                         final double[][] inverseOut,
                                                         final int dimensionUsed) {
        double logDet = 0.0;
        for (int i = 0; i < dimensionUsed; ++i) {
            for (int j = 0; j <= i; ++j) {
                double sum = cholesky[i][j];
                for (int k = 0; k < j; ++k) {
                    sum -= cholesky[i][k] * cholesky[j][k];
                }
                if (i == j) {
                    if (!(sum > 0.0) || !Double.isFinite(sum)) {
                        return Double.NaN;
                    }
                    final double diag = Math.sqrt(sum);
                    cholesky[i][j] = diag;
                    logDet += Math.log(diag);
                } else {
                    cholesky[i][j] = sum / cholesky[j][j];
                }
            }
            for (int j = i + 1; j < dimensionUsed; ++j) {
                cholesky[i][j] = 0.0;
            }
        }
        logDet *= 2.0;

        for (int column = 0; column < dimensionUsed; ++column) {
            for (int row = 0; row < dimensionUsed; ++row) {
                double sum = row == column ? 1.0 : 0.0;
                for (int k = 0; k < row; ++k) {
                    sum -= cholesky[row][k] * lowerInverse[k][column];
                }
                lowerInverse[row][column] = sum / cholesky[row][row];
            }
        }

        for (int i = 0; i < dimensionUsed; ++i) {
            for (int j = 0; j < dimensionUsed; ++j) {
                double sum = 0.0;
                for (int k = 0; k < dimensionUsed; ++k) {
                    sum += lowerInverse[k][i] * lowerInverse[k][j];
                }
                inverseOut[i][j] = sum;
            }
        }
        return logDet;
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

    private static double gershgorinLowerBound(final double[][] matrix, final int dimensionUsed) {
        double lowerBound = Double.POSITIVE_INFINITY;
        for (int i = 0; i < dimensionUsed; ++i) {
            double radius = 0.0;
            for (int j = 0; j < dimensionUsed; ++j) {
                if (i != j) {
                    radius += Math.abs(matrix[i][j]);
                }
            }
            lowerBound = Math.min(lowerBound, matrix[i][i] - radius);
        }
        return lowerBound;
    }

    private void emitSpdFailureDebugDump(final String context,
                                         final DenseMatrix64F originalSource,
                                         final DenseMatrix64F symmetrizedSource,
                                         final double jitterBase) {
        if (!Boolean.getBoolean(SPD_DEBUG_DUMP_PROPERTY)) {
            return;
        }
        final String path = "/tmp/ou_spd_failure_" + System.nanoTime() + ".json";
        try {
            final StringBuilder sb = new StringBuilder(8192);
            sb.append("{\n");
            sb.append("\"context\":\"").append(context).append("\",\n");
            sb.append("\"dimension\":").append(symmetrizedSource.numRows).append(",\n");
            sb.append("\"jitterBase\":").append(jitterBase).append(",\n");
            sb.append("\"minDiagonal\":").append(minDiagonal(symmetrizedSource)).append(",\n");
            sb.append("\"maxAbsDiagonal\":").append(maxAbsDiagonal(symmetrizedSource)).append(",\n");
            sb.append("\"gershgorinLowerBound\":").append(gershgorinLowerBound(symmetrizedSource)).append(",\n");
            sb.append("\"originalSource\":").append(jsonMatrix(originalSource)).append(",\n");
            sb.append("\"symmetrizedSource\":").append(jsonMatrix(symmetrizedSource)).append(",\n");
            sb.append("\"actualization\":").append(jsonMatrix(actualizationMatrix)).append(",\n");
            sb.append("\"displacement\":").append(jsonVector(displacementVector)).append(",\n");
            sb.append("\"branchPrecisionBuffer\":").append(jsonMatrix(branchPrecisionMatrix)).append(",\n");
            sb.append("\"branchCovarianceBuffer\":").append(jsonMatrix(branchCovarianceMatrix)).append(",\n");
            sb.append("\"aboveChildPrecisionBuffer\":").append(jsonMatrix(aboveChildPrecisionMatrix)).append(",\n");
            sb.append("\"aboveChildCovarianceBuffer\":").append(jsonMatrix(aboveChildCovarianceMatrix)).append(",\n");
            sb.append("\"aboveParentPrecisionBuffer\":").append(jsonMatrix(aboveParentPrecisionMatrix)).append(",\n");
            sb.append("\"aboveParentCovarianceBuffer\":").append(jsonMatrix(aboveParentCovarianceMatrix)).append("\n");
            sb.append("}\n");

            Files.write(Paths.get(path), sb.toString().getBytes(StandardCharsets.UTF_8));
            System.err.println("OU_SPD_FAILURE_DUMP " + path);
        } catch (final Exception e) {
            System.err.println("OU_SPD_FAILURE_DUMP_FAILED context=" + context + " reason=" + e.getMessage());
        }
    }

    private static String jsonMatrix(final DenseMatrix64F matrix) {
        final StringBuilder sb = new StringBuilder();
        sb.append("[");
        for (int i = 0; i < matrix.numRows; ++i) {
            if (i > 0) {
                sb.append(",");
            }
            sb.append("[");
            for (int j = 0; j < matrix.numCols; ++j) {
                if (j > 0) {
                    sb.append(",");
                }
                sb.append(matrix.unsafe_get(i, j));
            }
            sb.append("]");
        }
        sb.append("]");
        return sb.toString();
    }

    private static String jsonVector(final DenseMatrix64F vector) {
        final StringBuilder sb = new StringBuilder();
        sb.append("[");
        for (int i = 0; i < vector.numRows; ++i) {
            if (i > 0) {
                sb.append(",");
            }
            sb.append(vector.unsafe_get(i, 0));
        }
        sb.append("]");
        return sb.toString();
    }

    private static double minDiagonal(final DenseMatrix64F matrix) {
        final int d = Math.min(matrix.numRows, matrix.numCols);
        double min = Double.POSITIVE_INFINITY;
        for (int i = 0; i < d; ++i) {
            min = Math.min(min, matrix.unsafe_get(i, i));
        }
        return min;
    }

    private static double gershgorinLowerBound(final DenseMatrix64F matrix) {
        final int d = Math.min(matrix.numRows, matrix.numCols);
        double lowerBound = Double.POSITIVE_INFINITY;
        for (int i = 0; i < d; ++i) {
            double radius = 0.0;
            for (int j = 0; j < d; ++j) {
                if (i != j) {
                    radius += Math.abs(matrix.unsafe_get(i, j));
                }
            }
            lowerBound = Math.min(lowerBound, matrix.unsafe_get(i, i) - radius);
        }
        return lowerBound;
    }
}

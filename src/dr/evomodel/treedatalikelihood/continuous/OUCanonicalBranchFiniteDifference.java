package dr.evomodel.treedatalikelihood.continuous;

import dr.evomodel.continuous.ou.OUProcessModel;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalLocalTransitionAdjoints;
import dr.evomodel.treedatalikelihood.preorder.BranchSufficientStatistics;
import dr.inference.model.Parameter;
import org.ejml.data.DenseMatrix64F;

import java.util.Arrays;

final class OUCanonicalBranchFiniteDifference {

    private final OUGaussianBranchTransitionProvider branchTransitionProvider;
    private final OUCanonicalBranchWiring branchWiring;
    private final OUProcessModel processModel;
    private final CanonicalDebugOptions debugOptions;
    private final int dimension;
    private final OUCanonicalBranchGradientScratch scratch;

    OUCanonicalBranchFiniteDifference(final OUGaussianBranchTransitionProvider branchTransitionProvider,
                                      final OUCanonicalBranchWiring branchWiring,
                                      final OUProcessModel processModel,
                                      final CanonicalDebugOptions debugOptions,
                                      final int dimension,
                                      final OUCanonicalBranchGradientScratch scratch) {
        this.branchTransitionProvider = branchTransitionProvider;
        this.branchWiring = branchWiring;
        this.processModel = processModel;
        this.debugOptions = debugOptions;
        this.dimension = dimension;
        this.scratch = scratch;
    }

    double evaluateLocalSelectionScore(final double branchLength,
                                       final double[] optimum,
                                       final CanonicalLocalTransitionAdjoints adjoints) {
        refreshProcessSnapshots();
        return branchWiring.evaluateLocalSelectionScore(branchLength, optimum, adjoints);
    }

    double[] numericalLocalSelectionGradient(final double branchLength,
                                             final double[] optimum,
                                             final CanonicalLocalTransitionAdjoints adjoints,
                                             final Parameter requestedParameter,
                                             final RuntimeException cause) {
        final int parameterDimension = requestedParameter.getDimension();
        ensureParameterScratchCapacity(parameterDimension);
        for (int i = 0; i < parameterDimension; ++i) {
            scratch.parameterRestore[i] = requestedParameter.getParameterValue(i);
        }

        final double[] gradient = new double[parameterDimension];
        try {
            for (int i = 0; i < parameterDimension; ++i) {
                final double base = scratch.parameterRestore[i];
                final double step = 1.0e-6 * Math.max(1.0, Math.abs(base));

                requestedParameter.setParameterValue(i, base + step);
                final double plus = evaluateLocalSelectionScore(branchLength, optimum, adjoints);

                requestedParameter.setParameterValue(i, base - step);
                final double minus = evaluateLocalSelectionScore(branchLength, optimum, adjoints);

                if (!Double.isFinite(plus) || !Double.isFinite(minus)) {
                    throw new IllegalStateException(
                            "Non-finite local selection score in numerical fallback at index=" + i
                                    + ", base=" + base
                                    + ", step=" + step
                                    + ", plus=" + plus
                                    + ", minus=" + minus);
                }

                gradient[i] = (plus - minus) / (2.0 * step);
                if (!Double.isFinite(gradient[i])) {
                    throw new IllegalStateException(
                            "Non-finite local numerical gradient at index=" + i
                                    + ", plus=" + plus
                                    + ", minus=" + minus
                                    + ", step=" + step);
                }
                requestedParameter.setParameterValue(i, base);
            }
            return gradient;
        } catch (final RuntimeException e) {
            for (int i = 0; i < parameterDimension; ++i) {
                requestedParameter.setParameterValue(i, scratch.parameterRestore[i]);
            }
            final IllegalStateException wrapped = new IllegalStateException(
                    "Failed numerical local fallback after analytic native block selection gradient failure",
                    e);
            if (cause != null) {
                wrapped.addSuppressed(cause);
            }
            throw wrapped;
        } finally {
            for (int i = 0; i < parameterDimension; ++i) {
                if (requestedParameter.getParameterValue(i) != scratch.parameterRestore[i]) {
                    requestedParameter.setParameterValue(i, scratch.parameterRestore[i]);
                }
            }
        }
    }

    double[] numericalLocalSelectionGradientGeneric(final double branchLength,
                                                    final double[] optimum,
                                                    final CanonicalLocalTransitionAdjoints adjoints,
                                                    final Parameter requestedParameter) {
        final int parameterDimension = requestedParameter.getDimension();
        final double[] restore = new double[parameterDimension];
        for (int i = 0; i < parameterDimension; ++i) {
            restore[i] = requestedParameter.getParameterValue(i);
        }

        final double[] gradient = new double[parameterDimension];
        try {
            for (int i = 0; i < parameterDimension; ++i) {
                final double base = restore[i];
                final double step = 1.0e-6 * Math.max(1.0, Math.abs(base));

                requestedParameter.setParameterValue(i, base + step);
                final double plus = evaluateLocalSelectionScore(branchLength, optimum, adjoints);

                requestedParameter.setParameterValue(i, base - step);
                final double minus = evaluateLocalSelectionScore(branchLength, optimum, adjoints);

                gradient[i] = (plus - minus) / (2.0 * step);
                requestedParameter.setParameterValue(i, base);
            }
            return gradient;
        } finally {
            for (int i = 0; i < parameterDimension; ++i) {
                if (requestedParameter.getParameterValue(i) != restore[i]) {
                    requestedParameter.setParameterValue(i, restore[i]);
                }
            }
        }
    }

    double[] numericalLocalSelectionGradientGenericFromFrozenFactor(final double branchLength,
                                                                    final double[] optimum,
                                                                    final BranchSufficientStatistics statistics,
                                                                    final Parameter requestedParameter) {
        final int parameterDimension = requestedParameter.getDimension();
        final double[] restore = new double[parameterDimension];
        for (int i = 0; i < parameterDimension; ++i) {
            restore[i] = requestedParameter.getParameterValue(i);
        }

        final double[] gradient = new double[parameterDimension];
        try {
            for (int i = 0; i < parameterDimension; ++i) {
                final double base = restore[i];
                final double step = 1.0e-6 * Math.max(1.0, Math.abs(base));

                requestedParameter.setParameterValue(i, base + step);
                refreshProcessSnapshots();
                final double plus = branchWiring.evaluateFrozenLocalLogFactor(
                        branchLength, optimum, statistics.getAbove(), statistics.getBelow());

                requestedParameter.setParameterValue(i, base - step);
                refreshProcessSnapshots();
                final double minus = branchWiring.evaluateFrozenLocalLogFactor(
                        branchLength, optimum, statistics.getAbove(), statistics.getBelow());

                gradient[i] = (plus - minus) / (2.0 * step);
                requestedParameter.setParameterValue(i, base);
            }
            return gradient;
        } finally {
            for (int i = 0; i < parameterDimension; ++i) {
                if (requestedParameter.getParameterValue(i) != restore[i]) {
                    requestedParameter.setParameterValue(i, restore[i]);
                }
            }
            refreshProcessSnapshots();
        }
    }

    void emitDenseLocalSelectionDebug(final double branchLength,
                                      final double[] optimum,
                                      final CanonicalLocalTransitionAdjoints adjoints,
                                      final double[] analytic) {
        final double[] numeric = numericalLocalSelectionGradientDenseA(branchLength, optimum, adjoints);
        double maxAbs = 0.0;
        int maxIdx = -1;
        for (int i = 0; i < Math.min(analytic.length, numeric.length); ++i) {
            final double diff = Math.abs(analytic[i] - numeric[i]);
            if (diff > maxAbs) {
                maxAbs = diff;
                maxIdx = i;
            }
        }
        System.err.println("branchLocalSelectionFDDebug len=" + branchLength
                + " maxAbsDiff=" + maxAbs
                + " maxIdx=" + maxIdx
                + " analytic=" + Arrays.toString(analytic)
                + " numeric=" + Arrays.toString(numeric));
    }

    void maybeEmitLocalMeanFDDebug(final double branchLength,
                                   final double[] optimum,
                                   final BranchSufficientStatistics statistics,
                                   final int nodeNumber,
                                   final double[] analytic,
                                   final CanonicalLocalTransitionAdjoints localAdjoints) {
        if (!debugOptions.isBranchLocalMeanFiniteDifferenceEnabled()) {
            return;
        }
        final Integer requestedNode = debugOptions.getBranchLocalMeanFiniteDifferenceNode();
        if (requestedNode != null && nodeNumber != requestedNode) {
            return;
        }
        final double[] numericAdjointScore =
                numericalFrozenLocalMeanGradient(branchLength, optimum, localAdjoints);
        final double[] numericFrozenFactor =
                numericalFrozenLocalMeanGradientFromFrozenFactor(branchLength, optimum, statistics);
        final double baseAdjointScore = evaluateLocalSelectionScore(branchLength, optimum, localAdjoints);
        final double baseFrozenFactor = branchWiring.evaluateFrozenLocalLogFactor(
                branchLength,
                optimum,
                statistics.getAbove(),
                statistics.getAboveParent(),
                statistics.getBelow());
        double maxAbsAdjoint = 0.0;
        int maxIdxAdjoint = -1;
        double maxAbsFrozenFactor = 0.0;
        int maxIdxFrozenFactor = -1;
        for (int i = 0; i < analytic.length; ++i) {
            final double adjointDiff = Math.abs(analytic[i] - numericAdjointScore[i]);
            if (adjointDiff > maxAbsAdjoint) {
                maxAbsAdjoint = adjointDiff;
                maxIdxAdjoint = i;
            }
            final double frozenFactorDiff = Math.abs(analytic[i] - numericFrozenFactor[i]);
            if (frozenFactorDiff > maxAbsFrozenFactor) {
                maxAbsFrozenFactor = frozenFactorDiff;
                maxIdxFrozenFactor = i;
            }
        }
        System.err.println("branchLocalMeanFDDebug node=" + nodeNumber
                + " len=" + branchLength
                + " belowPrecisionHasInf=" + hasInfinity(statsBelowPrecision(statistics))
                + " aboveParentPrecisionHasInf=" + hasInfinity(statsAboveParentPrecision(statistics))
                + " belowHasMissingMask=" + hasMissingMask(statistics)
                + " maxAbsAdjointScoreDiff=" + maxAbsAdjoint
                + " maxIdxAdjointScore=" + maxIdxAdjoint
                + " maxAbsFrozenFactorDiff=" + maxAbsFrozenFactor
                + " maxIdxFrozenFactor=" + maxIdxFrozenFactor
                + " baseScoreDelta=" + (baseAdjointScore - baseFrozenFactor)
                + " analytic=" + Arrays.toString(analytic)
                + " numericAdjointScore=" + Arrays.toString(numericAdjointScore)
                + " numericFrozenFactor=" + Arrays.toString(numericFrozenFactor));
    }

    double[] numericalLocalSelectionGradientDenseA(final double branchLength,
                                                   final double[] optimum,
                                                   final CanonicalLocalTransitionAdjoints adjoints) {
        final int d = dimension;
        final double[] gradient = new double[d * d];
        final double[] restore = new double[d * d];
        final dr.inference.model.MatrixParameterInterface drift = processModel.getDriftMatrix();

        int k = 0;
        for (int row = 0; row < d; ++row) {
            for (int col = 0; col < d; ++col) {
                restore[k++] = drift.getParameterValue(row, col);
            }
        }

        try {
            k = 0;
            for (int row = 0; row < d; ++row) {
                for (int col = 0; col < d; ++col) {
                    final double base = restore[k];
                    final double step = 1.0e-6 * Math.max(1.0, Math.abs(base));
                    drift.setParameterValue(row, col, base + step);
                    final double plus = evaluateLocalSelectionScore(branchLength, optimum, adjoints);
                    drift.setParameterValue(row, col, base - step);
                    final double minus = evaluateLocalSelectionScore(branchLength, optimum, adjoints);
                    gradient[k] = (plus - minus) / (2.0 * step);
                    drift.setParameterValue(row, col, base);
                    ++k;
                }
            }
            return gradient;
        } finally {
            k = 0;
            for (int row = 0; row < d; ++row) {
                for (int col = 0; col < d; ++col) {
                    final double original = restore[k++];
                    if (drift.getParameterValue(row, col) != original) {
                        drift.setParameterValue(row, col, original);
                    }
                }
            }
        }
    }

    double[] numericalLocalSelectionGradientFromFrozenFactor(final double branchLength,
                                                             final double[] optimum,
                                                             final BranchSufficientStatistics statistics,
                                                             final CanonicalLocalTransitionAdjoints localAdjoints) {
        branchWiring.fillLocalAdjoints(branchLength, optimum, statistics, localAdjoints);
        return numericalLocalSelectionGradientDenseA(branchLength, optimum, localAdjoints);
    }

    private double[] numericalFrozenLocalMeanGradient(final double branchLength,
                                                      final double[] optimum,
                                                      final CanonicalLocalTransitionAdjoints localAdjoints) {
        final double[] gradient = new double[dimension];
        final double[] optimumWork = optimum.clone();
        for (int i = 0; i < dimension; ++i) {
            final double base = optimumWork[i];
            final double step = 1.0e-6 * Math.max(1.0, Math.abs(base));
            optimumWork[i] = base + step;
            final double plus = evaluateLocalSelectionScore(branchLength, optimumWork, localAdjoints);
            optimumWork[i] = base - step;
            final double minus = evaluateLocalSelectionScore(branchLength, optimumWork, localAdjoints);
            gradient[i] = (plus - minus) / (2.0 * step);
            optimumWork[i] = base;
        }
        return gradient;
    }

    private double[] numericalFrozenLocalMeanGradientFromFrozenFactor(final double branchLength,
                                                                      final double[] optimum,
                                                                      final BranchSufficientStatistics statistics) {
        final double[] gradient = new double[dimension];
        final double[] optimumWork = optimum.clone();
        for (int i = 0; i < dimension; ++i) {
            final double base = optimumWork[i];
            final double step = 1.0e-6 * Math.max(1.0, Math.abs(base));
            optimumWork[i] = base + step;
            final double plus = branchWiring.evaluateFrozenLocalLogFactor(
                    branchLength,
                    optimumWork,
                    statistics.getAbove(),
                    statistics.getAboveParent(),
                    statistics.getBelow());
            optimumWork[i] = base - step;
            final double minus = branchWiring.evaluateFrozenLocalLogFactor(
                    branchLength,
                    optimumWork,
                    statistics.getAbove(),
                    statistics.getAboveParent(),
                    statistics.getBelow());
            gradient[i] = (plus - minus) / (2.0 * step);
            optimumWork[i] = base;
        }
        return gradient;
    }

    private void ensureParameterScratchCapacity(final int parameterDimension) {
        if (scratch.parameterRestore.length < parameterDimension) {
            throw new IllegalStateException("parameter scratch too small for requested parameter dimension");
        }
    }

    private void refreshProcessSnapshots() {
        branchTransitionProvider.getProcessModel();
    }

    private static DenseMatrix64F statsBelowPrecision(final BranchSufficientStatistics statistics) {
        return statistics.getBelow() == null ? null : statistics.getBelow().getRawPrecision();
    }

    private static DenseMatrix64F statsAboveParentPrecision(final BranchSufficientStatistics statistics) {
        return statistics.getAboveParent() == null ? null : statistics.getAboveParent().getRawPrecision();
    }

    private static boolean hasMissingMask(final BranchSufficientStatistics statistics) {
        final int[] missing = statistics.getMissing();
        if (missing == null) {
            return false;
        }
        for (int value : missing) {
            if (value != 0) {
                return true;
            }
        }
        return false;
    }

    private static boolean hasInfinity(final DenseMatrix64F matrix) {
        if (matrix == null) {
            return false;
        }
        final double[] data = matrix.getData();
        for (double value : data) {
            if (Double.isInfinite(value)) {
                return true;
            }
        }
        return false;
    }
}

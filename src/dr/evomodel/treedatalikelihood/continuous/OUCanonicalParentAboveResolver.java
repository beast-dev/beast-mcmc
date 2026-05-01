package dr.evomodel.treedatalikelihood.continuous;

import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalNumerics;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalNumericsOptions;
import dr.evomodel.treedatalikelihood.preorder.NormalSufficientStatistics;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

final class OUCanonicalParentAboveResolver {

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

    private final int dimension;
    private final DenseMatrix64F actualizationMatrix;
    private final DenseMatrix64F displacementVector;
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
    private final double[][] symmetric2DScratch;
    private final double[][] adjusted2DScratch;
    private final double[][] choleskyScratch;
    private final double[][] lowerInverseScratch;
    private final CanonicalNumericsOptions numericsOptions;
    private final OUCanonicalBranchSpdDebugContext spdDebugContext;
    private int parentAboveCompareReportCount = 0;

    OUCanonicalParentAboveResolver(final int dimension,
                                   final DenseMatrix64F actualizationMatrix,
                                   final DenseMatrix64F displacementVector,
                                   final DenseMatrix64F branchPrecisionMatrix,
                                   final DenseMatrix64F branchCovarianceMatrix,
                                   final CanonicalNumericsOptions numericsOptions) {
        this.dimension = dimension;
        this.actualizationMatrix = actualizationMatrix;
        this.displacementVector = displacementVector;
        this.branchCovarianceMatrix = branchCovarianceMatrix;
        this.aboveChildPrecisionMatrix = new DenseMatrix64F(dimension, dimension);
        this.aboveChildCovarianceMatrix = new DenseMatrix64F(dimension, dimension);
        this.aboveParentPrecisionMatrix = new DenseMatrix64F(dimension, dimension);
        this.aboveParentCovarianceMatrix = new DenseMatrix64F(dimension, dimension);
        this.actualizationInverseMatrix = new DenseMatrix64F(dimension, dimension);
        this.aboveParentMeanVector = new DenseMatrix64F(dimension, 1);
        this.centeredChildMeanVector = new DenseMatrix64F(dimension, 1);
        this.scratchMatrix = new DenseMatrix64F(dimension, dimension);
        this.secondaryScratchMatrix = new DenseMatrix64F(dimension, dimension);
        this.symmetric2DScratch = new double[dimension][dimension];
        this.adjusted2DScratch = new double[dimension][dimension];
        this.choleskyScratch = new double[dimension][dimension];
        this.lowerInverseScratch = new double[dimension][dimension];
        this.numericsOptions = numericsOptions;
        this.spdDebugContext = new OUCanonicalBranchSpdDebugContext(
                actualizationMatrix,
                displacementVector,
                branchPrecisionMatrix,
                branchCovarianceMatrix,
                aboveChildPrecisionMatrix,
                aboveChildCovarianceMatrix,
                aboveParentPrecisionMatrix,
                aboveParentCovarianceMatrix);
    }

    DenseMatrix64F recoverOrUseParentAbovePrecision(final NormalSufficientStatistics aboveChild,
                                                    final NormalSufficientStatistics aboveParent) {
        if (Boolean.getBoolean(FORCE_RECOVER_PARENT_ABOVE_FROM_CHILD_PROPERTY)) {
            return recoverParentAbovePrecision(aboveChild);
        }
        if (aboveParent == null) {
            throw new IllegalStateException(
                    "Missing parent-above message for canonical branch contribution."
                            + " childAbovePrecisionFinite=" + CanonicalNumerics.isFinite(aboveChild.getRawPrecision())
                            + " childAboveMeanFinite=" + CanonicalNumerics.isFinite(aboveChild.getRawMean()));
        }
        if (Boolean.getBoolean(COMPARE_PARENT_ABOVE_RECOVERED_PROPERTY)) {
            compareParentAboveAgainstRecovered(aboveChild, aboveParent);
        }
        final DenseMatrix64F parentPrecisionRaw = aboveParent.getRawPrecision();
        final DenseMatrix64F parentMeanRaw = aboveParent.getRawMean();
        final boolean parentMeanFinite = CanonicalNumerics.isFinite(parentMeanRaw);
        final boolean parentPrecisionHasNaN = CanonicalNumerics.hasNaN(parentPrecisionRaw);
        final boolean parentPrecisionHasInfinity = CanonicalNumerics.hasInfinity(parentPrecisionRaw);
        if (!parentMeanFinite || parentPrecisionHasNaN) {
            throw new IllegalStateException(
                    "Invalid parent-above message."
                            + " parentMeanFinite=" + parentMeanFinite
                            + " parentPrecisionHasNaN=" + parentPrecisionHasNaN
                            + " parentPrecisionHasInfinity=" + parentPrecisionHasInfinity
                            + " parentPrecisionSummary=" + CanonicalNumerics.summarizeDenseMatrix(parentPrecisionRaw)
                            + " parentMeanSummary=" + CanonicalNumerics.summarizeDenseMatrix(parentMeanRaw)
                            + " childAbovePrecisionFinite=" + CanonicalNumerics.isFinite(aboveChild.getRawPrecision())
                            + " childAboveMeanFinite=" + CanonicalNumerics.isFinite(aboveChild.getRawMean())
                            + " childAbovePrecisionSummary=" + CanonicalNumerics.summarizeDenseMatrix(aboveChild.getRawPrecision())
                            + " childAboveMeanSummary=" + CanonicalNumerics.summarizeDenseMatrix(aboveChild.getRawMean()));
        }
        if (parentPrecisionHasInfinity) {
            if (!Boolean.parseBoolean(System.getProperty(ALLOW_INFINITE_PARENT_ABOVE_RECOVERY_PROPERTY, "true"))) {
                throw new IllegalStateException(
                        "Parent-above precision contains infinite entries and infinite-recovery is disabled."
                                + " parentPrecisionSummary=" + CanonicalNumerics.summarizeDenseMatrix(parentPrecisionRaw)
                                + " parentMeanSummary=" + CanonicalNumerics.summarizeDenseMatrix(parentMeanRaw));
            }
            // Exact parent constraints need a finite precision for local canonical algebra.
            // Keep the supplied parent mean fixed; only the precision is recovered.
            final DenseMatrix64F recoveredPrecision = recoverParentAbovePrecision(aboveChild);
            aboveParentMeanVector.set(parentMeanRaw);
            return recoveredPrecision;
        }
        aboveParentPrecisionMatrix.set(parentPrecisionRaw);
        aboveParentMeanVector.set(parentMeanRaw);
        return aboveParentPrecisionMatrix;
    }

    DenseMatrix64F requireParentAbovePrecision(final NormalSufficientStatistics aboveParent) {
        if (aboveParent == null) {
            throw new IllegalStateException("Missing parent-above message for canonical branch contribution.");
        }
        final DenseMatrix64F parentPrecisionRaw = aboveParent.getRawPrecision();
        final DenseMatrix64F parentMeanRaw = aboveParent.getRawMean();
        final boolean parentMeanFinite = CanonicalNumerics.isFinite(parentMeanRaw);
        final boolean parentPrecisionHasNaN = CanonicalNumerics.hasNaN(parentPrecisionRaw);
        final boolean parentPrecisionHasInfinity = CanonicalNumerics.hasInfinity(parentPrecisionRaw);
        if (!parentMeanFinite || parentPrecisionHasNaN || parentPrecisionHasInfinity) {
            throw new IllegalStateException(
                    "Invalid parent-above message."
                            + " parentMeanFinite=" + parentMeanFinite
                            + " parentPrecisionHasNaN=" + parentPrecisionHasNaN
                            + " parentPrecisionHasInfinity=" + parentPrecisionHasInfinity
                            + " parentPrecisionSummary=" + CanonicalNumerics.summarizeDenseMatrix(parentPrecisionRaw)
                            + " parentMeanSummary=" + CanonicalNumerics.summarizeDenseMatrix(parentMeanRaw));
        }
        aboveParentPrecisionMatrix.set(parentPrecisionRaw);
        aboveParentMeanVector.set(parentMeanRaw);
        return aboveParentPrecisionMatrix;
    }

    DenseMatrix64F recoverParentAbovePrecision(final NormalSufficientStatistics above) {
        aboveChildPrecisionMatrix.set(above.getRawPrecision());
        safeInvertSymmetricPositiveDefinite(
                aboveChildPrecisionMatrix,
                aboveChildCovarianceMatrix,
                "recoverParentAbovePrecision:aboveChildPrecision->aboveChildCovariance");

        for (int i = 0; i < dimension; ++i) {
            centeredChildMeanVector.unsafe_set(i, 0,
                    above.getRawMean().unsafe_get(i, 0) - displacementVector.unsafe_get(i, 0));
        }

        secondaryScratchMatrix.set(actualizationMatrix);
        final double solveJitter = numericsOptions.jitterBase(
                CanonicalNumerics.maxAbsDiagonal(secondaryScratchMatrix));
        CanonicalNumerics.addDiagonalJitter(secondaryScratchMatrix, solveJitter);
        safeInvert(secondaryScratchMatrix, actualizationInverseMatrix);
        CommonOps.mult(actualizationInverseMatrix, centeredChildMeanVector, aboveParentMeanVector);

        aboveParentCovarianceMatrix.set(aboveChildCovarianceMatrix);
        CommonOps.subtractEquals(aboveParentCovarianceMatrix, branchCovarianceMatrix);
        CanonicalNumerics.symmetrizeInPlace(aboveParentCovarianceMatrix);
        CommonOps.mult(actualizationInverseMatrix, aboveParentCovarianceMatrix, scratchMatrix);
        CommonOps.multTransB(scratchMatrix, actualizationInverseMatrix, aboveParentCovarianceMatrix);
        CanonicalNumerics.symmetrizeInPlace(aboveParentCovarianceMatrix);
        safeInvertSymmetricPositiveDefinite(
                aboveParentCovarianceMatrix,
                aboveParentPrecisionMatrix,
                "recoverParentAbovePrecision:aboveParentCovariance->aboveParentPrecision");
        return aboveParentPrecisionMatrix;
    }

    DenseMatrix64F getAboveParentMeanVector() {
        return aboveParentMeanVector;
    }

    CanonicalNumerics.DenseSpdFailureDump getSpdFailureDump() {
        return spdDebugContext;
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

    private void safeInvertSymmetricPositiveDefinite(final DenseMatrix64F source,
                                                     final DenseMatrix64F inverseOut,
                                                     final String context) {
        CanonicalNumerics.safeInvertSymmetricPositiveDefinite(
                source,
                inverseOut,
                scratchMatrix,
                secondaryScratchMatrix,
                symmetric2DScratch,
                adjusted2DScratch,
                choleskyScratch,
                lowerInverseScratch,
                numericsOptions,
                context,
                spdDebugContext);
    }

    private void safeInvert(final DenseMatrix64F source,
                            final DenseMatrix64F inverseOut) {
        CanonicalNumerics.safeInvertWithJitter(
                source, inverseOut, secondaryScratchMatrix, numericsOptions);
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
}

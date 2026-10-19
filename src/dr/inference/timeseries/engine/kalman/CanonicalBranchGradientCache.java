package dr.inference.timeseries.engine.kalman;

import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalBranchMessageContribution;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalLocalTransitionAdjoints;
import dr.evomodel.treedatalikelihood.continuous.canonical.math.MatrixOps;
import dr.inference.timeseries.core.TimeGrid;
import dr.inference.timeseries.representation.CachedGaussianTransitionRepresentation;
import dr.inference.timeseries.representation.GaussianTransitionRepresentation;
import dr.inference.timeseries.representation.TransitionMomentsView;

import java.util.Arrays;

/**
 * Per-smoother cache of branch-local canonical adjoints shared by gradient formulas.
 */
public final class CanonicalBranchGradientCache {

    private final CanonicalBranchMessageContribution[] contributions;
    private final boolean[] contributionKnown;
    private final CanonicalLocalTransitionAdjoints[] adjoints;
    private final double[] transitionMatrices;
    private final MatrixSlice[] transitionMatrixViews;
    private final GaussianTransitionRepresentation transitionRepresentation;
    private final CachedGaussianTransitionRepresentation cachedTransitionRepresentation;
    private final TimeGrid timeGrid;
    private final MomentForwardTrajectory momentTrajectory;
    private final MomentBranchSmootherStats[] smootherStats;
    private final TransitionMomentsView transitionMomentsView;
    private final CanonicalBranchMessageContribution workingContribution;
    private final double[] transitionMatrixWorkspace;
    private final double[] covarianceDifference;
    private final double[] matrixScratch1;
    private final double[] matrixScratch2;
    private final double[] meanDifference;
    private SharedCanonicalTimeSeriesSchedule sharedSchedule;
    private CanonicalForwardTrajectory knownTrajectory;
    private boolean known;
    private long buildCount;
    private final int matrixLength;

    CanonicalBranchGradientCache(final int timeCount,
                                 final int stateDimension,
                                 final GaussianTransitionRepresentation transitionRepresentation,
                                 final TimeGrid timeGrid,
                                 final MomentForwardTrajectory momentTrajectory,
                                 final MomentBranchSmootherStats[] smootherStats) {
        final int branchCount = Math.max(0, timeCount - 1);
        this.matrixLength = stateDimension * stateDimension;
        this.contributions = new CanonicalBranchMessageContribution[branchCount];
        this.contributionKnown = new boolean[branchCount];
        this.adjoints = new CanonicalLocalTransitionAdjoints[branchCount];
        this.transitionMatrices = new double[branchCount * matrixLength];
        this.transitionMatrixViews = new MatrixSlice[branchCount];
        this.transitionRepresentation = transitionRepresentation;
        this.cachedTransitionRepresentation =
                transitionRepresentation instanceof CachedGaussianTransitionRepresentation
                        ? (CachedGaussianTransitionRepresentation) transitionRepresentation
                        : null;
        this.timeGrid = timeGrid;
        this.momentTrajectory = momentTrajectory;
        this.smootherStats = smootherStats;
        this.transitionMomentsView = new TransitionMomentsView();
        this.workingContribution = new CanonicalBranchMessageContribution(stateDimension);
        for (int i = 0; i < branchCount; ++i) {
            adjoints[i] = new CanonicalLocalTransitionAdjoints(stateDimension);
            transitionMatrixViews[i] = new MatrixSlice();
        }
        this.transitionMatrixWorkspace = new double[matrixLength];
        this.covarianceDifference = new double[matrixLength];
        this.matrixScratch1 = new double[matrixLength];
        this.matrixScratch2 = new double[matrixLength];
        this.meanDifference = new double[stateDimension];
        this.sharedSchedule = null;
        this.knownTrajectory = null;
        this.known = false;
        this.buildCount = 0L;
    }

    public void setSharedSchedule(final SharedCanonicalTimeSeriesSchedule sharedSchedule) {
        this.sharedSchedule = sharedSchedule;
        makeDirty();
    }

    public void makeDirty() {
        known = false;
        knownTrajectory = null;
        Arrays.fill(contributionKnown, false);
    }

    public void ensure(final CanonicalForwardTrajectory trajectory) {
        if (known) {
            return;
        }
        ++buildCount;
        knownTrajectory = trajectory;
        Arrays.fill(contributionKnown, false);
        for (int t = 0; t < adjoints.length; ++t) {
            final int transitionMatrixOffset = transitionMatrixOffset(t);
            if (sharedSchedule != null) {
                fillContributionFromRtsMoments(t, workingContribution);
            }
            fillAdjointsFromRtsMoments(t, adjoints[t]);
            if (cachedTransitionRepresentation != null
                    && cachedTransitionRepresentation.getTransitionMomentsView(
                    t,
                    t + 1,
                    timeGrid,
                    transitionMomentsView)) {
                transitionMatrixViews[t].set(transitionMomentsView.getTransitionMatrix(), 0);
            } else {
                transitionRepresentation.getTransitionMatrixFlat(
                        t, t + 1, timeGrid, transitionMatrixWorkspace);
                System.arraycopy(transitionMatrixWorkspace, 0,
                        transitionMatrices, transitionMatrixOffset, matrixLength);
                setOwnedTransitionMatrix(t, transitionMatrixOffset);
            }
        }
        known = true;
    }

    public CanonicalBranchMessageContribution getContribution(final int branchIndex) {
        if (!known || knownTrajectory == null) {
            throw new IllegalStateException("Canonical branch gradient cache is not initialized");
        }
        if (contributions[branchIndex] == null) {
            contributions[branchIndex] =
                    new CanonicalBranchMessageContribution(adjoints[branchIndex].getDimension());
        }
        if (!contributionKnown[branchIndex]) {
            fillContributionFromRtsMoments(branchIndex, contributions[branchIndex]);
            contributionKnown[branchIndex] = true;
        }
        return contributions[branchIndex];
    }

    public CanonicalLocalTransitionAdjoints getAdjoints(final int branchIndex) {
        return adjoints[branchIndex];
    }

    public double[] getTransitionMatrixData(final int branchIndex) {
        return transitionMatrixViews[branchIndex].data;
    }

    public int getTransitionMatrixOffset(final int branchIndex) {
        return transitionMatrixViews[branchIndex].offset;
    }

    public long getBuildCount() {
        return buildCount;
    }

    private void fillAdjointsFromRtsMoments(final int branchIndex,
                                            final CanonicalLocalTransitionAdjoints out) {
        final int d = out.getDimension();
        final int nextTime = branchIndex + 1;
        final int nextVectorOffset = momentTrajectory.stateVectorOffset(nextTime);
        final int nextMatrixOffset = momentTrajectory.stateMatrixOffset(nextTime);
        final double[] predictedPrecision = knownTrajectory.predictedStates[nextTime].precision;

        for (int i = 0; i < d; ++i) {
            meanDifference[i] = smootherStats[nextTime].smoothedMean[i]
                    - momentTrajectory.predictedMeans[nextVectorOffset + i];
        }
        MatrixOps.matVec(predictedPrecision, meanDifference, out.dLogL_df, d);

        for (int i = 0; i < matrixLength; ++i) {
            covarianceDifference[i] = smootherStats[nextTime].smoothedCovariance[i]
                    - momentTrajectory.predictedCovariances[nextMatrixOffset + i];
        }
        MatrixOps.matMulTransposedRight(
                covarianceDifference,
                smootherStats[branchIndex].smootherGain,
                matrixScratch1,
                d);
        for (int i = 0; i < d; ++i) {
            final int rowOffset = i * d;
            for (int j = 0; j < d; ++j) {
                matrixScratch1[rowOffset + j] +=
                        meanDifference[i] * smootherStats[branchIndex].smoothedMean[j];
            }
        }
        MatrixOps.matMul(predictedPrecision, matrixScratch1, out.dLogL_dF, d);

        System.arraycopy(
                smootherStats[nextTime].smoothedCovariance,
                0,
                matrixScratch1,
                0,
                matrixLength);
        for (int i = 0; i < d; ++i) {
            final int rowOffset = i * d;
            for (int j = 0; j < d; ++j) {
                matrixScratch1[rowOffset + j] += meanDifference[i] * meanDifference[j];
            }
        }
        MatrixOps.matMul(predictedPrecision, matrixScratch1, matrixScratch2, d);
        MatrixOps.matMul(matrixScratch2, predictedPrecision, out.dLogL_dOmega, d);
        for (int i = 0; i < matrixLength; ++i) {
            out.dLogL_dOmega[i] = 0.5
                    * (out.dLogL_dOmega[i] - predictedPrecision[i]);
        }
        MatrixOps.symmetrize(out.dLogL_dOmega, d);
    }

    private void fillContributionFromRtsMoments(final int branchIndex,
                                                final CanonicalBranchMessageContribution out) {
        final int d = out.getDimension();
        final MomentBranchSmootherStats current = smootherStats[branchIndex];
        final MomentBranchSmootherStats next = smootherStats[branchIndex + 1];
        if (sharedSchedule != null
                && sharedSchedule.fillContributionFromRtsMoments(
                branchIndex, current, next, out)) {
            return;
        }
        MatrixOps.matMul(
                current.smootherGain,
                next.smoothedCovariance,
                matrixScratch1,
                d);
        for (int i = 0; i < d; ++i) {
            final double currentMeanI = current.smoothedMean[i];
            final double nextMeanI = next.smoothedMean[i];
            out.dLogL_dInformationX[i] = currentMeanI;
            out.dLogL_dInformationY[i] = nextMeanI;
            for (int j = 0; j < d; ++j) {
                final int ij = i * d + j;
                out.dLogL_dPrecisionXX[ij] = -0.5
                        * (current.smoothedCovariance[ij]
                        + currentMeanI * current.smoothedMean[j]);
                out.dLogL_dPrecisionXY[ij] = -0.5
                        * (matrixScratch1[ij]
                        + currentMeanI * next.smoothedMean[j]);
                out.dLogL_dPrecisionYX[ij] = -0.5
                        * (matrixScratch1[j * d + i]
                        + nextMeanI * current.smoothedMean[j]);
                out.dLogL_dPrecisionYY[ij] = -0.5
                        * (next.smoothedCovariance[ij]
                        + nextMeanI * next.smoothedMean[j]);
            }
        }
        out.dLogL_dLogNormalizer = -1.0;
    }

    private int transitionMatrixOffset(final int branchIndex) {
        return branchIndex * matrixLength;
    }

    private void setOwnedTransitionMatrix(final int branchIndex,
                                          final int transitionMatrixOffset) {
        transitionMatrixViews[branchIndex].set(transitionMatrices, transitionMatrixOffset);
    }

    private static final class MatrixSlice {
        private double[] data;
        private int offset;

        private void set(final double[] data,
                         final int offset) {
            this.data = data;
            this.offset = offset;
        }
    }
}

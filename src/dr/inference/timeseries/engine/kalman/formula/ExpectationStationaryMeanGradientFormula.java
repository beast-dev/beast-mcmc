package dr.inference.timeseries.engine.kalman.formula;

import dr.inference.model.MatrixParameter;
import dr.inference.model.Parameter;
import dr.inference.timeseries.core.TimeGrid;
import dr.evomodel.continuous.ou.OUProcessModel;
import dr.evomodel.continuous.ou.blockdiagonal.BlockDiagonalCanonicalParameterization;
import dr.evomodel.treedatalikelihood.continuous.canonical.math.MatrixOps;
import dr.inference.timeseries.engine.kalman.BranchSmootherStats;
import dr.inference.timeseries.engine.kalman.ForwardTrajectory;
import dr.inference.timeseries.representation.GaussianTransitionRepresentation;

/**
 * Analytical gradient of the log-likelihood with respect to the OU stationary/optimal mean.
 *
 * <p>This parameter affects:
 * <ul>
 *   <li>the branch offset {@code f_t = (I - F_t) mu}</li>
 *   <li>the initial latent mean {@code x_0 ~ N(mu, P0)}</li>
 * </ul>
 * so the total gradient is the sum of branch contributions and the initial-state term.
 */
public final class ExpectationStationaryMeanGradientFormula implements ExpectationGradientFormula {

    private final Parameter stationaryMeanParameter;
    private final MatrixParameter initialCovarianceParameter;
    private final int stateDimension;
    private final OUProcessModel processModel;

    private final ExpectationGaussianBranchGradientAdjoints branchAdjoints;
    private final double[] initialCovarianceFlat;
    private final double[] initialCovInv;
    private final double[] initialCovCholesky;
    private final double[] initialCovLowerInverse;
    private final double[] initialMeanAdjoint;
    private final double[] currentMean;
    private final double[] stateDiff;
    private final double[] denseGradient;

    public ExpectationStationaryMeanGradientFormula(final Parameter stationaryMeanParameter,
                                         final MatrixParameter initialCovarianceParameter,
                                         final int stateDimension) {
        this(null, stationaryMeanParameter, initialCovarianceParameter, stateDimension);
    }

    public ExpectationStationaryMeanGradientFormula(final OUProcessModel processModel,
                                         final Parameter stationaryMeanParameter,
                                         final MatrixParameter initialCovarianceParameter,
                                         final int stateDimension) {
        if (stationaryMeanParameter == null) {
            throw new IllegalArgumentException("stationaryMeanParameter must not be null");
        }
        if (initialCovarianceParameter == null) {
            throw new IllegalArgumentException("initialCovarianceParameter must not be null");
        }
        if (stateDimension < 1) {
            throw new IllegalArgumentException("stateDimension must be at least 1");
        }
        this.processModel = processModel;
        this.stationaryMeanParameter = stationaryMeanParameter;
        this.initialCovarianceParameter = initialCovarianceParameter;
        this.stateDimension = stateDimension;

        this.branchAdjoints = new ExpectationGaussianBranchGradientAdjoints(stateDimension);
        final int matrixSize = stateDimension * stateDimension;
        this.initialCovarianceFlat = new double[matrixSize];
        this.initialCovInv = new double[matrixSize];
        this.initialCovCholesky = new double[matrixSize];
        this.initialCovLowerInverse = new double[matrixSize];
        this.initialMeanAdjoint = new double[stateDimension];
        this.currentMean = new double[stateDimension];
        this.stateDiff = new double[stateDimension];
        this.denseGradient = new double[stateDimension];
    }

    @Override
    public boolean supportsParameter(final Parameter parameter) {
        return parameter == stationaryMeanParameter;
    }

    @Override
    public double[] computeGradient(final Parameter parameter,
                                    final BranchSmootherStats[] smootherStats,
                                    final ForwardTrajectory trajectory,
                                    final GaussianTransitionRepresentation repr,
                                    final TimeGrid timeGrid) {
        if (stationaryMeanParameter.getDimension() == 1) {
            return computeScalarGradient(smootherStats, trajectory, timeGrid);
        }

        BlockDiagonalFormulaSupport.zero(denseGradient);
        BlockDiagonalFormulaSupport.fillCurrentMean(
                stationaryMeanParameter,
                stateDimension,
                currentMean);

        final int timeCount = trajectory.timeCount;
        final BlockDiagonalCanonicalParameterization blockParameterization =
                BlockDiagonalFormulaSupport.canonicalParameterization(processModel);
        for (int t = 0; t < timeCount - 1; ++t) {
            final BranchSmootherStats curr = smootherStats[t];
            final BranchSmootherStats next = smootherStats[t + 1];
            final int transitionMatrixOffset = trajectory.branchMatrixOffset(t);
            final int transitionOffsetOffset = trajectory.branchVectorOffset(t);

            branchAdjoints.computeMeanAdjoint(
                    curr,
                    next,
                    trajectory.transitionMatrices,
                    transitionMatrixOffset,
                    trajectory.transitionOffsets,
                    transitionOffsetOffset,
                    trajectory.stepCovariances,
                    transitionMatrixOffset);

            if (blockParameterization != null) {
                blockParameterization.accumulateMeanGradient(
                        timeGrid.getDelta(t, t + 1),
                        branchAdjoints.dLogL_df(),
                        denseGradient);
            } else {
                BlockDiagonalFormulaSupport.accumulateBranchMeanGradientFlat(
                        trajectory.transitionMatrices,
                        transitionMatrixOffset,
                        branchAdjoints.dLogL_df(),
                        denseGradient);
            }
        }

        fillInitialCovarianceFlat();
        if (!MatrixOps.tryCholesky(initialCovarianceFlat, initialCovCholesky, stateDimension)) {
            throw new IllegalArgumentException("Matrix is not positive definite");
        }
        MatrixOps.invertFromCholesky(
                initialCovCholesky, initialCovLowerInverse, initialCovInv, stateDimension);
        for (int i = 0; i < stateDimension; ++i) {
            stateDiff[i] = smootherStats[0].smoothedMean[i] - currentMean[i];
        }
        MatrixOps.matVec(initialCovInv, stateDiff, initialMeanAdjoint, stateDimension);
        for (int i = 0; i < stateDimension; ++i) {
            denseGradient[i] += initialMeanAdjoint[i];
        }

        return BlockDiagonalFormulaSupport.projectMeanGradient(
                stationaryMeanParameter,
                stateDimension,
                denseGradient);
    }

    private double[] computeScalarGradient(final BranchSmootherStats[] smootherStats,
                                           final ForwardTrajectory trajectory,
                                           final TimeGrid timeGrid) {
        double scalarGradient = 0.0;
        final double meanValue = stationaryMeanParameter.getParameterValue(0);
        final int timeCount = trajectory.timeCount;
        final BlockDiagonalCanonicalParameterization blockParameterization =
                BlockDiagonalFormulaSupport.canonicalParameterization(processModel);

        for (int t = 0; t < timeCount - 1; ++t) {
            final BranchSmootherStats curr = smootherStats[t];
            final BranchSmootherStats next = smootherStats[t + 1];
            final int transitionMatrixOffset = trajectory.branchMatrixOffset(t);
            final int transitionOffsetOffset = trajectory.branchVectorOffset(t);

            branchAdjoints.computeMeanAdjoint(
                    curr,
                    next,
                    trajectory.transitionMatrices,
                    transitionMatrixOffset,
                    trajectory.transitionOffsets,
                    transitionOffsetOffset,
                    trajectory.stepCovariances,
                    transitionMatrixOffset);

            if (blockParameterization != null) {
                scalarGradient += blockParameterization.accumulateScalarMeanGradient(
                        timeGrid.getDelta(t, t + 1),
                        branchAdjoints.dLogL_df());
            } else {
                scalarGradient += BlockDiagonalFormulaSupport.accumulateScalarBranchMeanGradientFlat(
                        trajectory.transitionMatrices,
                        transitionMatrixOffset,
                        branchAdjoints.dLogL_df());
            }
        }

        fillInitialCovarianceFlat();
        if (!MatrixOps.tryCholesky(initialCovarianceFlat, initialCovCholesky, stateDimension)) {
            throw new IllegalArgumentException("Matrix is not positive definite");
        }
        MatrixOps.invertFromCholesky(
                initialCovCholesky, initialCovLowerInverse, initialCovInv, stateDimension);
        for (int i = 0; i < stateDimension; ++i) {
            stateDiff[i] = smootherStats[0].smoothedMean[i] - meanValue;
        }
        MatrixOps.matVec(initialCovInv, stateDiff, initialMeanAdjoint, stateDimension);
        for (int i = 0; i < stateDimension; ++i) {
            scalarGradient += initialMeanAdjoint[i];
        }

        return new double[]{scalarGradient};
    }

    private void fillInitialCovarianceFlat() {
        for (int i = 0; i < stateDimension; ++i) {
            final int rowOffset = i * stateDimension;
            for (int j = 0; j < stateDimension; ++j) {
                initialCovarianceFlat[rowOffset + j] =
                        initialCovarianceParameter.getParameterValue(i, j);
            }
        }
    }
}

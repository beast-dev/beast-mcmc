package dr.evomodel.continuous.ou.diagonal;

import dr.evomodel.continuous.MultivariateElasticModel;
import dr.evomodel.continuous.ou.kernel.DiagonalOUBranchTransformKernel;
import dr.evomodel.continuous.ou.kernel.OUBranchKernelWorkspace;
import dr.evomodel.continuous.ou.kernel.OUPreparedBranch;
import dr.evomodel.continuous.ou.math.OUPositiveDefiniteMatrixInversionOptions;
import dr.evomodel.continuous.ou.math.PositiveDefiniteMatrixOps;
import dr.inference.model.DiagonalMatrix;
import dr.inference.model.MatrixParameterInterface;

public final class DiagonalOUBranchKernel implements DiagonalOUBranchTransformKernel {

    private final DiagonalMatrix selectionMatrix;
    private final int dimension;

    public DiagonalOUBranchKernel(final DiagonalMatrix selectionMatrix) {
        if (selectionMatrix == null) {
            throw new IllegalArgumentException("selectionMatrix must not be null");
        }
        this.selectionMatrix = selectionMatrix;
        this.dimension = selectionMatrix.getRowDimension();
    }

    @Override
    public int getDimension() {
        return dimension;
    }

    @Override
    public DiagonalOUPreparedBranch createPreparedBranch() {
        return new DiagonalOUPreparedBranch(dimension);
    }

    @Override
    public DiagonalOUBranchKernelWorkspace createWorkspace() {
        return new DiagonalOUBranchKernelWorkspace(dimension);
    }

    @Override
    public void prepareBranch(final double effectiveBranchLength,
                              final MatrixParameterInterface diffusionCovariance,
                              final double[] stationaryMean,
                              final OUPreparedBranch preparedBranch,
                              final OUBranchKernelWorkspace workspace) {
        if (!(preparedBranch instanceof DiagonalOUPreparedBranch)) {
            throw new IllegalArgumentException("preparedBranch must be DiagonalOUPreparedBranch");
        }
        if (!(workspace instanceof DiagonalOUBranchKernelWorkspace)) {
            throw new IllegalArgumentException("workspace must be DiagonalOUBranchKernelWorkspace");
        }
        if (stationaryMean == null || stationaryMean.length < dimension) {
            throw new IllegalArgumentException("stationaryMean length must be at least " + dimension);
        }
        final DiagonalOUPreparedBranch prepared = (DiagonalOUPreparedBranch) preparedBranch;
        final DiagonalOUBranchKernelWorkspace work = (DiagonalOUBranchKernelWorkspace) workspace;

        fillDiffusionCovariance(diffusionCovariance, work.branchCovarianceScratch);
        for (int i = 0; i < dimension; ++i) {
            final double alpha = selectionMatrix.getParameterValue(i, i);
            final double scaled = -alpha * effectiveBranchLength;
            prepared.actualizationDiagonal()[i] = Math.exp(scaled);
            work.oneMinusActualization[i] = -Math.expm1(scaled);
            prepared.displacement()[i] = work.oneMinusActualization[i] * stationaryMean[i];
        }

        for (int i = 0; i < dimension; ++i) {
            final double alphaI = selectionMatrix.getParameterValue(i, i);
            final int iOff = i * dimension;
            for (int j = 0; j < dimension; ++j) {
                final double alphaJ = selectionMatrix.getParameterValue(j, j);
                final double denominator = alphaI + alphaJ;
                final double qij = work.branchCovarianceScratch[iOff + j];
                final double covariance;
                if (denominator == 0.0) {
                    covariance = effectiveBranchLength * qij;
                } else {
                    final double factor = work.oneMinusActualization[i]
                            + work.oneMinusActualization[j]
                            - work.oneMinusActualization[i] * work.oneMinusActualization[j];
                    covariance = qij * factor / denominator;
                }
                work.branchCovariance[iOff + j] = covariance;
            }
        }

        final double logDet = PositiveDefiniteMatrixOps.copyAndInvertPositiveDefiniteFlat(
                work.branchCovariance,
                dimension,
                work.branchCovarianceScratch,
                prepared.precision(),
                work.branchCholeskyScratch,
                work.branchLowerInverseScratch,
                OUPositiveDefiniteMatrixInversionOptions.OU_TREE);
        final double logNormalizer = 0.5 * (dimension * MultivariateElasticModel.LOG2PI
                + logDet
                + quadratic(prepared.precision(), prepared.displacement()));
        prepared.markCurrent(effectiveBranchLength, logNormalizer);
    }

    @Override
    public double[] actualizationDiagonal(final OUPreparedBranch preparedBranch) {
        return asDiagonalPreparedBranch(preparedBranch).actualizationDiagonal();
    }

    @Override
    public double[] displacement(final OUPreparedBranch preparedBranch) {
        return asDiagonalPreparedBranch(preparedBranch).displacement();
    }

    @Override
    public double[] precision(final OUPreparedBranch preparedBranch) {
        return asDiagonalPreparedBranch(preparedBranch).precision();
    }

    @Override
    public double logNormalizer(final OUPreparedBranch preparedBranch) {
        return asDiagonalPreparedBranch(preparedBranch).logNormalizer();
    }

    @Override
    public boolean isCurrent(final OUPreparedBranch preparedBranch, final double branchLength) {
        return asDiagonalPreparedBranch(preparedBranch).isCurrent(branchLength);
    }

    private void fillDiffusionCovariance(final MatrixParameterInterface diffusionMatrix,
                                         final double[] out) {
        for (int i = 0; i < dimension; ++i) {
            final int iOff = i * dimension;
            for (int j = 0; j < dimension; ++j) {
                out[iOff + j] = diffusionMatrix.getParameterValue(i, j);
            }
        }
    }

    private double quadratic(final double[] precision, final double[] value) {
        double out = 0.0;
        for (int i = 0; i < dimension; ++i) {
            double row = 0.0;
            final int iOff = i * dimension;
            for (int j = 0; j < dimension; ++j) {
                row += precision[iOff + j] * value[j];
            }
            out += value[i] * row;
        }
        return out;
    }

    private DiagonalOUPreparedBranch asDiagonalPreparedBranch(final OUPreparedBranch preparedBranch) {
        if (!(preparedBranch instanceof DiagonalOUPreparedBranch)) {
            throw new IllegalArgumentException("preparedBranch must be DiagonalOUPreparedBranch");
        }
        return (DiagonalOUPreparedBranch) preparedBranch;
    }
}

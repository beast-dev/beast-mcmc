package dr.inference.timeseries.engine.gaussian;

import dr.evomodel.continuous.ou.OUProcessModel;
import dr.evomodel.continuous.ou.orthogonalblockdiagonal.BlockDiagonalCanonicalParameterization;
import dr.evomodel.continuous.ou.orthogonalblockdiagonal.BlockDiagonalNativeCanonicalParameterization;
import dr.evomodel.treedatalikelihood.continuous.canonical.gradient.CanonicalSelectionGradientProjector;
import dr.inference.model.AbstractBlockDiagonalTwoByTwoMatrixParameter;
import dr.inference.model.OrthogonalMatrixProvider;
import dr.inference.model.Parameter;
import dr.inference.timeseries.representation.GaussianTransitionRepresentation;
import dr.inference.timeseries.representation.KernelBackedGaussianTransitionRepresentation;

/**
 * Shared block-diagonal plumbing used by the dense and canonical gradient formulas.
 */
final class BlockDiagonalFormulaSupport {

    private static final String DISABLE_ORTHOGONAL_NATIVE_SELECTION_PATH_PROPERTY =
            "beast.experimental.disableOrthogonalNativeSelectionPath";
    private static final String DISABLE_BLOCK_DIAGONAL_NATIVE_SELECTION_PATH_PROPERTY =
            "beast.experimental.disableBlockDiagonalNativeSelectionPath";

    private BlockDiagonalFormulaSupport() { }

    static boolean supportsSelectionParameter(final Parameter selectionMatrixParameter,
                                              final Parameter requestedParameter) {
        if (requestedParameter == selectionMatrixParameter) {
            return true;
        }
        if (!(selectionMatrixParameter instanceof AbstractBlockDiagonalTwoByTwoMatrixParameter)) {
            return false;
        }
        return supportsBlockSelectionSubParameter(
                (AbstractBlockDiagonalTwoByTwoMatrixParameter) selectionMatrixParameter,
                requestedParameter);
    }

    static boolean supportsNativeSelectionParameter(final AbstractBlockDiagonalTwoByTwoMatrixParameter blockParameter,
                                                    final Parameter requestedParameter) {
        if (requestedParameter == blockParameter) {
            return false;
        }
        return supportsBlockSelectionSubParameter(blockParameter, requestedParameter);
    }

    static boolean isNativeSelectionAvailable(final OUProcessModel processModel,
                                              final Parameter selectionMatrixParameter) {
        return processModel != null
                && selectionMatrixParameter instanceof AbstractBlockDiagonalTwoByTwoMatrixParameter
                && processModel.getCovarianceGradientMethod()
                == OUProcessModel.CovarianceGradientMethod.STATIONARY_LYAPUNOV
                && processModel.getSelectionMatrixParameterization()
                instanceof BlockDiagonalNativeCanonicalParameterization;
    }

    static boolean shouldUseNativeSelectionPath(final Parameter selectionMatrixParameter,
                                                final Parameter requestedParameter,
                                                final OUProcessModel processModel) {
        return !nativeSelectionPathDisabled()
                && requestedParameter != selectionMatrixParameter
                && isNativeSelectionAvailable(processModel, selectionMatrixParameter);
    }

    static OUProcessModel ouProcessModel(final GaussianTransitionRepresentation repr) {
        if (repr instanceof OUProcessModel) {
            return (OUProcessModel) repr;
        }
        if (repr instanceof KernelBackedGaussianTransitionRepresentation) {
            return ((KernelBackedGaussianTransitionRepresentation) repr).getProcessModel();
        }
        return null;
    }

    static BlockDiagonalNativeCanonicalParameterization nativeParameterization(
            final OUProcessModel processModel) {
        if (processModel != null
                && processModel.getSelectionMatrixParameterization()
                instanceof BlockDiagonalNativeCanonicalParameterization) {
            return (BlockDiagonalNativeCanonicalParameterization)
                    processModel.getSelectionMatrixParameterization();
        }
        return null;
    }

    static BlockDiagonalCanonicalParameterization canonicalParameterization(
            final OUProcessModel processModel) {
        if (processModel != null
                && processModel.getSelectionMatrixParameterization()
                instanceof BlockDiagonalCanonicalParameterization) {
            return (BlockDiagonalCanonicalParameterization)
                    processModel.getSelectionMatrixParameterization();
        }
        return null;
    }

    static double[] assembleNativeSelectionGradient(
            final int stateDimension,
            final Parameter requestedParameter,
            final AbstractBlockDiagonalTwoByTwoMatrixParameter blockParameter,
            final double[] nativeGradient,
            final double[] gradientR) {
        return CanonicalSelectionGradientProjector.assembleBlockGradientResultFlat(
                stateDimension,
                requestedParameter,
                blockParameter,
                nativeGradient,
                gradientR);
    }

    static double[] pullBackDenseSelectionGradient(
            final int stateDimension,
            final Parameter requestedParameter,
            final AbstractBlockDiagonalTwoByTwoMatrixParameter blockParameter,
            final double[] denseGradient,
            final CanonicalSelectionGradientProjector.Workspace workspace) {
        return CanonicalSelectionGradientProjector.pullBackDenseGradientToBlock(
                stateDimension,
                requestedParameter,
                blockParameter,
                denseGradient,
                workspace);
    }

    static void fillCurrentMean(final Parameter stationaryMeanParameter,
                                final int stateDimension,
                                final double[] out) {
        if (stationaryMeanParameter.getDimension() == 1) {
            final double value = stationaryMeanParameter.getParameterValue(0);
            for (int i = 0; i < stateDimension; ++i) {
                out[i] = value;
            }
        } else if (stationaryMeanParameter.getDimension() == stateDimension) {
            for (int i = 0; i < stateDimension; ++i) {
                out[i] = stationaryMeanParameter.getParameterValue(i);
            }
        } else {
            throw new IllegalStateException("stationaryMean dimension must be 1 or stateDimension");
        }
    }

    static double[] projectMeanGradient(final Parameter stationaryMeanParameter,
                                        final int stateDimension,
                                        final double[] denseSource) {
        if (stationaryMeanParameter.getDimension() == stateDimension) {
            return denseSource.clone();
        }
        double sum = 0.0;
        for (int i = 0; i < stateDimension; ++i) {
            sum += denseSource[i];
        }
        return new double[]{sum};
    }

    static void accumulateBranchMeanGradient(final double[][] transitionMatrix,
                                             final double[] dLogL_df,
                                             final double[] accumulator) {
        final int d = dLogL_df.length;
        for (int j = 0; j < d; ++j) {
            double sum = dLogL_df[j];
            for (int i = 0; i < d; ++i) {
                sum -= transitionMatrix[i][j] * dLogL_df[i];
            }
            accumulator[j] += sum;
        }
    }

    static void accumulateBranchMeanGradientFlat(final double[] transitionMatrix,
                                                 final double[] dLogL_df,
                                                 final double[] accumulator) {
        final int d = dLogL_df.length;
        for (int j = 0; j < d; ++j) {
            double sum = dLogL_df[j];
            for (int i = 0; i < d; ++i) {
                sum -= transitionMatrix[i * d + j] * dLogL_df[i];
            }
            accumulator[j] += sum;
        }
    }

    static double accumulateScalarBranchMeanGradient(final double[][] transitionMatrix,
                                                     final double[] dLogL_df) {
        final int d = dLogL_df.length;
        double accumulator = 0.0;
        for (int j = 0; j < d; ++j) {
            double sum = dLogL_df[j];
            for (int i = 0; i < d; ++i) {
                sum -= transitionMatrix[i][j] * dLogL_df[i];
            }
            accumulator += sum;
        }
        return accumulator;
    }

    static double accumulateScalarBranchMeanGradientFlat(final double[] transitionMatrix,
                                                         final double[] dLogL_df) {
        final int d = dLogL_df.length;
        double accumulator = 0.0;
        for (int j = 0; j < d; ++j) {
            double sum = dLogL_df[j];
            for (int i = 0; i < d; ++i) {
                sum -= transitionMatrix[i * d + j] * dLogL_df[i];
            }
            accumulator += sum;
        }
        return accumulator;
    }

    static void zero(final double[] vector) {
        for (int i = 0; i < vector.length; ++i) {
            vector[i] = 0.0;
        }
    }

    private static boolean nativeSelectionPathDisabled() {
        return Boolean.getBoolean(DISABLE_BLOCK_DIAGONAL_NATIVE_SELECTION_PATH_PROPERTY)
                || Boolean.getBoolean(DISABLE_ORTHOGONAL_NATIVE_SELECTION_PATH_PROPERTY);
    }

    private static boolean supportsBlockSelectionSubParameter(
            final AbstractBlockDiagonalTwoByTwoMatrixParameter blockParameter,
            final Parameter requestedParameter) {
        if (requestedParameter == blockParameter.getParameter()) {
            return true;
        }
        if (requestedParameter == blockParameter.getRotationMatrixParameter()) {
            return true;
        }
        if (blockParameter.getRotationMatrixParameter() instanceof OrthogonalMatrixProvider
                && requestedParameter == ((OrthogonalMatrixProvider) blockParameter.getRotationMatrixParameter())
                .getOrthogonalParameter()) {
            return true;
        }
        if (requestedParameter == blockParameter.getScalarBlockParameter()
                && blockParameter.getScalarBlockParameter().getDimension() > 0) {
            return true;
        }
        for (int i = 0; i < blockParameter.getTwoByTwoParameterFamilyCount(); ++i) {
            if (requestedParameter == blockParameter.getTwoByTwoBlockParameter(i)) {
                return true;
            }
        }
        return false;
    }
}

package dr.evomodel.treedatalikelihood.continuous.canonical.adapter;

import dr.evomodel.continuous.ou.kernel.DiagonalOUBranchTransformKernel;
import dr.evomodel.continuous.ou.kernel.OUBranchKernelWorkspace;
import dr.evomodel.continuous.ou.kernel.OUPreparedBranch;
import dr.evomodel.treedatalikelihood.continuous.ContinuousBranchTimeModel;
import dr.evomodel.treedatalikelihood.continuous.canonical.CanonicalNativeBranchMessageProvider;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianMessageOps;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianState;
import dr.evomodel.treedatalikelihood.continuous.canonical.workspace.BranchGradientWorkspace;
import dr.evomodel.treedatalikelihood.continuous.observationmodel.TipObservationMode;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Parameter;

final class DiagonalOUCanonicalBranchMessageAdapter implements CanonicalNativeBranchMessageProvider {

    private final DiagonalOUBranchTransformKernel kernel;
    private final int dimension;
    private final MatrixParameterInterface diffusionMatrix;
    private final Parameter stationaryMean;
    private final ContinuousBranchTimeModel branchTimeModel;
    private final OUPreparedBranch[] cache;
    private final ThreadLocal<Workspace> workspaceLocal;

    DiagonalOUCanonicalBranchMessageAdapter(final DiagonalOUBranchTransformKernel kernel,
                                            final int nodeCount,
                                            final MatrixParameterInterface diffusionMatrix,
                                            final Parameter stationaryMean,
                                            final ContinuousBranchTimeModel branchTimeModel) {
        if (kernel == null) {
            throw new IllegalArgumentException("kernel must not be null");
        }
        if (diffusionMatrix == null) {
            throw new IllegalArgumentException("diffusionMatrix must not be null");
        }
        if (stationaryMean == null) {
            throw new IllegalArgumentException("stationaryMean must not be null");
        }
        if (branchTimeModel == null) {
            throw new IllegalArgumentException("branchTimeModel must not be null");
        }
        this.kernel = kernel;
        this.dimension = kernel.getDimension();
        this.diffusionMatrix = diffusionMatrix;
        this.stationaryMean = stationaryMean;
        this.branchTimeModel = branchTimeModel;
        this.cache = new OUPreparedBranch[nodeCount];
        this.workspaceLocal = ThreadLocal.withInitial(() -> new Workspace(kernel));
    }

    @Override
    public boolean fillBackwardMessageForState(final int childNodeIndex,
                                               final CanonicalGaussianState childState,
                                               final BranchGradientWorkspace branchWorkspace,
                                               final CanonicalGaussianState out) {
        final OUPreparedBranch entry = entry(childNodeIndex);
        if (entry == null) {
            return false;
        }
        CanonicalGaussianMessageOps.pushBackwardDiagonalTransition(
                childState,
                kernel.actualizationDiagonal(entry),
                kernel.displacement(entry),
                kernel.precision(entry),
                kernel.logNormalizer(entry),
                branchWorkspace.gaussianWorkspace,
                out);
        return true;
    }

    @Override
    public boolean fillBackwardMessageForObservedTip(final int childNodeIndex,
                                                     final double[] observedValue,
                                                     final CanonicalGaussianState out) {
        final OUPreparedBranch entry = entry(childNodeIndex);
        if (entry == null) {
            return false;
        }
        CanonicalGaussianMessageOps.conditionOnObservedSecondBlockDiagonalTransition(
                kernel.actualizationDiagonal(entry),
                kernel.displacement(entry),
                kernel.precision(entry),
                kernel.logNormalizer(entry),
                observedValue,
                out);
        return true;
    }

    @Override
    public boolean fillForwardMessage(final int childNodeIndex,
                                      final CanonicalGaussianState parentState,
                                      final BranchGradientWorkspace branchWorkspace,
                                      final CanonicalGaussianState out) {
        final OUPreparedBranch entry = entry(childNodeIndex);
        if (entry == null) {
            return false;
        }
        CanonicalGaussianMessageOps.pushForwardDiagonalTransition(
                parentState,
                kernel.actualizationDiagonal(entry),
                kernel.displacement(entry),
                kernel.precision(entry),
                kernel.logNormalizer(entry),
                branchWorkspace.gaussianWorkspace,
                out);
        return true;
    }

    @Override
    public boolean hasNativeBranchMessages() {
        return true;
    }

    @Override
    public boolean supportsPostorderInternalMessage() {
        return true;
    }

    @Override
    public boolean supportsPostorderTipMessage(final TipObservationMode mode) {
        return mode == TipObservationMode.EXACT_IDENTITY
                || mode == TipObservationMode.MISSING;
    }

    @Override
    public boolean supportsPreorderForwardMessage() {
        return true;
    }

    @Override
    public void clearBranchMessageCache() {
        for (OUPreparedBranch entry : cache) {
            if (entry != null) {
                entry.markInvalid();
            }
        }
    }

    private OUPreparedBranch entry(final int childNodeIndex) {
        final double branchLength = branchTimeModel.getEffectiveBranchTime(childNodeIndex);
        if (!(branchLength > 0.0) || !Double.isFinite(branchLength)) {
            return null;
        }

        OUPreparedBranch entry = cache[childNodeIndex];
        if (entry == null) {
            entry = kernel.createPreparedBranch();
            cache[childNodeIndex] = entry;
        }
        if (!kernel.isCurrent(entry, branchLength)) {
            final Workspace workspace = workspaceLocal.get();
            fillStationaryMean(workspace.stationaryMean);
            kernel.prepareBranch(
                    branchLength,
                    diffusionMatrix,
                    workspace.stationaryMean,
                    entry,
                    workspace.kernelWorkspace);
        }
        return entry;
    }

    private void fillStationaryMean(final double[] out) {
        if (stationaryMean.getDimension() == 1) {
            final double value = stationaryMean.getParameterValue(0);
            for (int i = 0; i < dimension; ++i) {
                out[i] = value;
            }
            return;
        }
        for (int i = 0; i < dimension; ++i) {
            out[i] = stationaryMean.getParameterValue(i);
        }
    }

    private static final class Workspace {
        private final OUBranchKernelWorkspace kernelWorkspace;
        private final double[] stationaryMean;

        private Workspace(final DiagonalOUBranchTransformKernel kernel) {
            this.kernelWorkspace = kernel.createWorkspace();
            this.stationaryMean = new double[kernel.getDimension()];
        }
    }
}

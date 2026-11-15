package dr.evomodel.treedatalikelihood.continuousNew.canonical.adapter;

import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.continuous.MultivariateDiffusionModel;
import dr.evomodel.continuous.MultivariateElasticModel;
import dr.evomodel.continuous.ou.OUProcessModel;
import dr.evomodel.continuous.ou.blockdiagonal.BlockBasisCanonicalTransitionBuilder;
import dr.evomodel.continuous.ou.blockdiagonal.BlockBasisCanonicalAdjointBuilder;
import dr.evomodel.continuous.ou.blockdiagonal.BlockBasisGlobalSnapshot;
import dr.evomodel.treedatalikelihood.continuous.ContinuousRateTransformation;
import dr.evomodel.treedatalikelihood.continuousNew.CanonicalOUOptions;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.CanonicalOUTransitionProvider;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.CanonicalPreparedBranchSnapshot;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.CanonicalPreparedBranchSnapshotProvider;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.CanonicalTransitionCacheDiagnostics;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.CanonicalTransitionPreloader;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.gradient.CanonicalSelectionGradientPullback;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.gradient.HomogeneousBlockBasisGradientCapability;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.message.CanonicalGaussianTransition;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.message.CanonicalLocalTransitionAdjoints;
import dr.inference.model.AbstractBlockDiagonalTwoByTwoMatrixParameter;
import dr.inference.model.AbstractModel;
import dr.inference.model.MatrixParameter;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.util.TaskPool;
import org.ejml.data.DenseMatrix64F;

import java.util.Arrays;

/**
 * Homogeneous OU transition provider that keeps all branch transitions in the
 * block basis {@code z = R^{-1}(y - beta)}.
 */
final class HomogeneousBlockBasisCanonicalOUBranchTransitionProvider extends AbstractModel
        implements CanonicalOUTransitionProvider,
        CanonicalPreparedBranchSnapshotProvider,
        CanonicalTransitionCacheDiagnostics,
        CanonicalTransitionPreloader,
        HomogeneousBlockBasisGradientCapability {

    private final Tree tree;
    private final MultivariateDiffusionModel diffusionModel;
    private final MultivariateElasticModel elasticModel;
    private final AbstractBlockDiagonalTwoByTwoMatrixParameter blockParameter;
    private final MatrixParameter diffusionCovariance;
    private final MatrixParameter blockBasisDiffusionCovariance;
    private final CanonicalOUDiffusionSnapshot diffusionSnapshot;
    private final CanonicalBranchLengthSnapshot branchLengths;
    private final BlockBasisGlobalSnapshot snapshot;
    private final BlockBasisCanonicalTransitionBuilder transitionBuilder;
    private final BlockBasisCanonicalAdjointBuilder adjointBuilder;
    private final OUProcessModel processModel;
    private final CanonicalTransitionCacheEntry[] entries;
    private final double[][] transitionMomentMatrices;
    private final double[][] transitionMomentCovariances;
    private final CanonicalTransitionCacheEpoch epoch;
    private final CanonicalTransitionCacheDiagnosticsRecorder diagnostics;
    private final double[] branchLengthTransitionDerivative;
    private final double[] branchLengthCovarianceDerivative;

    private boolean dirty = true;
    private CanonicalTransitionCacheInvalidationReason dirtyReason =
            CanonicalTransitionCacheInvalidationReason.MODEL_CHANGED;

    static boolean isSupported(final MultivariateElasticModel elasticModel) {
        return elasticModel.getStrengthOfSelectionMatrixParameter()
                instanceof AbstractBlockDiagonalTwoByTwoMatrixParameter;
    }

    HomogeneousBlockBasisCanonicalOUBranchTransitionProvider(final Tree tree,
                                                            final MultivariateElasticModel elasticModel,
                                                            final MultivariateDiffusionModel diffusionModel,
                                                            final BranchRateModel rateModel,
                                                            final ContinuousRateTransformation rateTransformation) {
        super("homogeneousBlockBasisCanonicalOUTransitionProvider");
        if (!isSupported(elasticModel)) {
            throw new UnsupportedOperationException(
                    "Homogeneous block-basis canonical OU requires a block-diagonal selection parametrization.");
        }
        this.tree = tree;
        this.elasticModel = elasticModel;
        this.diffusionModel = diffusionModel;
        this.blockParameter = (AbstractBlockDiagonalTwoByTwoMatrixParameter)
                elasticModel.getStrengthOfSelectionMatrixParameter();
        this.diffusionCovariance = new MatrixParameter(
                "blockBasisCanonicalOuProvider.diffusion",
                diffusionModel.getDimension(),
                diffusionModel.getDimension());
        this.blockBasisDiffusionCovariance = new MatrixParameter(
                "blockBasisCanonicalOuProvider.diffusionBlockBasis",
                diffusionModel.getDimension(),
                diffusionModel.getDimension());
        this.diffusionSnapshot = new CanonicalOUDiffusionSnapshot(
                diffusionModel.getDimension(),
                diffusionModel,
                diffusionCovariance);
        this.branchLengths = new CanonicalBranchLengthSnapshot(tree, rateModel, rateTransformation);
        this.snapshot = new BlockBasisGlobalSnapshot(blockParameter, diffusionCovariance);
        this.transitionBuilder = new BlockBasisCanonicalTransitionBuilder(snapshot.getStructure());
        this.adjointBuilder = new BlockBasisCanonicalAdjointBuilder(snapshot.getStructure());
        this.entries = new CanonicalTransitionCacheEntry[tree.getNodeCount()];
        this.transitionMomentMatrices = new double[tree.getNodeCount()][];
        this.transitionMomentCovariances = new double[tree.getNodeCount()][];
        for (int i = 0; i < entries.length; ++i) {
            entries[i] = new CanonicalTransitionCacheEntry();
            entries[i].processIndex = 0;
        }
        this.epoch = new CanonicalTransitionCacheEpoch();
        this.diagnostics = new CanonicalTransitionCacheDiagnosticsRecorder(
                CanonicalTransitionCacheOptions.fromSystemProperties().isDiagnosticsEnabled());
        this.branchLengthTransitionDerivative =
                new double[diffusionModel.getDimension() * diffusionModel.getDimension()];
        this.branchLengthCovarianceDerivative =
                new double[diffusionModel.getDimension() * diffusionModel.getDimension()];

        final MatrixParameter initialCovariance = new MatrixParameter(
                "blockBasisCanonicalOuProvider.initial",
                diffusionModel.getDimension(),
                diffusionModel.getDimension());
        setIdentity(initialCovariance);
        final Parameter zeroMean = new Parameter.Default(diffusionModel.getDimension());
        this.processModel = new OUProcessModel(
                "blockBasisCanonicalOuProvider.process",
                diffusionModel.getDimension(),
                blockParameter,
                blockBasisDiffusionCovariance,
                zeroMean,
                initialCovariance);

        addModel(diffusionModel);
        addModel(elasticModel);
        addVariable(blockParameter);
    }

    @Override
    public CanonicalSelectionGradientPullback createBlockBasisSelectionGradientPullback(final int dimension) {
        return new BlockBasisCanonicalSelectionGradientPullback(this, dimension);
    }

    /** Rotation {@code R} of the current block basis, {@code y = R z + beta}. */
    DenseMatrix64F getBlockBasisRotation() {
        return getBlockBasisSnapshot().getRotation();
    }

    /** Inverse rotation {@code R^{-1}} of the current block basis. */
    DenseMatrix64F getBlockBasisRotationInverse() {
        return getBlockBasisSnapshot().getRotationInverse();
    }

    /** {@code log |det R|} of the current block basis. */
    double getBlockBasisLogAbsDetRotation() {
        return getBlockBasisSnapshot().getLogAbsDetR();
    }

    public BlockBasisGlobalSnapshot getBlockBasisSnapshot() {
        ensureCurrentSnapshot();
        return snapshot;
    }

    public BlockBasisCanonicalAdjointBuilder getBlockBasisAdjointBuilder() {
        return adjointBuilder;
    }

    public AbstractBlockDiagonalTwoByTwoMatrixParameter getBlockParameter() {
        return blockParameter;
    }

    @Override
    public double contractBranchLengthGradient(final double branchLength,
                                               final CanonicalLocalTransitionAdjoints localAdjoints) {
        ensureCurrentSnapshot();
        synchronized (this) {
            transitionBuilder.fillBranchLengthDerivatives(
                    snapshot,
                    branchLength,
                    branchLengthTransitionDerivative,
                    branchLengthCovarianceDerivative);
            double score = 0.0;
            for (int i = 0; i < branchLengthTransitionDerivative.length; ++i) {
                score += localAdjoints.dLogL_dF[i] * branchLengthTransitionDerivative[i];
                score += localAdjoints.dLogL_dOmega[i] * branchLengthCovarianceDerivative[i];
            }
            return score;
        }
    }

    void ensureCurrentSnapshot() {
        if (!dirty) {
            return;
        }
        diagnostics.recordSnapshotRefresh();
        diffusionSnapshot.refresh();
        snapshot.ensureCurrent();
        copyBlockBasisDiffusionCovariance();
        clearEntries(dirtyReason);
        dirty = false;
        dirtyReason = CanonicalTransitionCacheInvalidationReason.MODEL_CHANGED;
    }

    @Override
    public int getDimension() {
        return diffusionModel.getDimension();
    }

    @Override
    public void fillCanonicalTransition(final int childNodeIndex,
                                        final CanonicalGaussianTransition out) {
        out.copyFrom(ensureTransition(childNodeIndex));
    }

    @Override
    public CanonicalGaussianTransition getCanonicalTransitionView(final int childNodeIndex) {
        return ensureTransition(childNodeIndex);
    }

    @Override
    public void preloadCanonicalTransitions(final int rootNodeIndex,
                                            final TaskPool taskPool) {
        if (taskPool == null || taskPool.getNumThreads() <= 1) {
            return;
        }
        taskPool.fork(
                (childNodeIndex, thread) -> {
                    if (childNodeIndex != rootNodeIndex) {
                        ensureTransition(childNodeIndex);
                    }
                });
    }

    @Override
    public double getEffectiveBranchLength(final int childNodeIndex) {
        return branchLengths.getEffectiveBranchLength(childNodeIndex);
    }

    @Override
    public void fillTraitCovarianceFlat(final double[] out) {
        ensureCurrentSnapshot();
        System.arraycopy(snapshot.getQDBasis().data, 0, out, 0, out.length);
    }

    @Override
    public OUProcessModel getProcessModel() {
        ensureCurrentSnapshot();
        return processModel;
    }

    @Override
    public void fillCanonicalTransitionForLength(final double branchLength,
                                                 final CanonicalGaussianTransition out) {
        ensureCurrentSnapshot();
        synchronized (this) {
            transitionBuilder.fillCanonicalTransition(snapshot, branchLength, out);
        }
    }

    @Override
    public CanonicalPreparedBranchSnapshot getPreparedBranchSnapshot(final int childNodeIndex) {
        ensureTransition(childNodeIndex);
        return entries[childNodeIndex].snapshot;
    }

    @Override
    public void fillTransitionMatrixFlat(final double branchLength, final double[] out) {
        ensureCurrentSnapshot();
        synchronized (this) {
            transitionBuilder.fillTransitionMatrix(snapshot, branchLength, out);
        }
    }

    @Override
    public void fillTransitionOffset(final double branchLength, final double[] out) {
        Arrays.fill(out, 0.0);
    }

    @Override
    public void fillTransitionCovarianceFlat(final double branchLength, final double[] out) {
        ensureCurrentSnapshot();
        synchronized (this) {
            transitionBuilder.fillTransitionCovariance(snapshot, branchLength, out);
        }
    }

    @Override
    public void fillTransitionMomentsFlat(final int childNodeIndex,
                                          final double branchLength,
                                          final double[] transitionMatrixOut,
                                          final double[] transitionOffsetOut,
                                          final double[] transitionCovarianceOut) {
        final CanonicalTransitionCacheEntry entry = ensureTransitionEntry(childNodeIndex);
        if (Double.doubleToLongBits(entry.effectiveBranchLength) != Double.doubleToLongBits(branchLength)) {
            throw new IllegalStateException("Cached block-basis moments for child " + childNodeIndex
                    + " have effective length " + entry.effectiveBranchLength
                    + " but caller requested " + branchLength + ".");
        }
        System.arraycopy(
                transitionMomentMatrices[childNodeIndex],
                0,
                transitionMatrixOut,
                0,
                getDimension() * getDimension());
        System.arraycopy(
                transitionMomentCovariances[childNodeIndex],
                0,
                transitionCovarianceOut,
                0,
                getDimension() * getDimension());
        Arrays.fill(transitionOffsetOut, 0.0);
    }

    @Override
    public boolean useTransitionMomentsForBranchAdjoints(final int childNodeIndex) {
        return true;
    }

    @Override
    public String pushDiagnosticPhase(final String phase) {
        ensureCurrentSnapshot();
        branchLengths.beginPhase();
        return diagnostics.pushPhase(phase);
    }

    @Override
    public void popDiagnosticPhase(final String previousPhase) {
        diagnostics.popPhase(previousPhase);
        branchLengths.endPhase();
    }

    @Override
    public void reportTransitionCacheDiagnostics(final String label) {
        diagnostics.report(label);
    }

    @Override
    public long getTransitionCacheMissCount(final String phase) {
        return diagnostics.phaseMisses(phase);
    }

    @Override
    public long getTransitionCacheRequestCount(final String phase) {
        return diagnostics.phaseRequests(phase);
    }

    @Override
    public long getTransitionCacheRebuildCount() {
        return diagnostics.transitionRebuilds();
    }

    @Override
    public long getPreparedBranchRebuildCount() {
        return diagnostics.preparedBasisRebuilds();
    }

    @Override
    public long getTransitionCacheClearCount(final String reason) {
        return diagnostics.clearCount(CanonicalTransitionCacheInvalidationReason.valueOf(reason));
    }

    @Override
    protected void handleModelChangedEvent(final Model model, final Object object, final int index) {
        if (model == diffusionModel) {
            markDirty(CanonicalTransitionCacheInvalidationReason.DIFFUSION_CHANGED);
        } else if (model == elasticModel) {
            markDirty(CanonicalTransitionCacheInvalidationReason.SELECTION_CHANGED);
        } else {
            markDirty(CanonicalTransitionCacheInvalidationReason.MODEL_CHANGED);
        }
    }

    @Override
    protected void handleVariableChangedEvent(final Variable variable,
                                              final int index,
                                              final Parameter.ChangeType type) {
        if (variable == blockParameter) {
            markDirty(CanonicalTransitionCacheInvalidationReason.SELECTION_CHANGED);
        } else {
            markDirty(CanonicalTransitionCacheInvalidationReason.MODEL_CHANGED);
        }
    }

    @Override
    public void storeState() {
        diagnostics.recordStore();
    }

    @Override
    public void restoreState() {
        diagnostics.recordRestore();
        clearEntries(CanonicalTransitionCacheInvalidationReason.RESTORE_STATE);
        dirty = true;
        dirtyReason = CanonicalTransitionCacheInvalidationReason.RESTORE_STATE;
    }

    @Override
    public void acceptState() {
        diagnostics.recordAccept();
    }

    private CanonicalGaussianTransition ensureTransition(final int childNodeIndex) {
        return ensureTransitionEntry(childNodeIndex).transition;
    }

    private CanonicalTransitionCacheEntry ensureTransitionEntry(final int childNodeIndex) {
        ensureCurrentSnapshot();
        diagnostics.recordRequest(diagnostics.currentPhase());
        final CanonicalTransitionCacheEntry entry = entries[childNodeIndex];
        final double effectiveBranchLength = getEffectiveBranchLength(childNodeIndex);
        if (!entry.isTransitionCurrent(effectiveBranchLength, epoch)) {
            synchronized (this) {
                if (!entry.isTransitionCurrent(effectiveBranchLength, epoch)) {
                    diagnostics.recordMiss(diagnostics.currentPhase());
                    diagnostics.recordTransitionRebuild();
                    if (entry.hasBranchLengthChanged(effectiveBranchLength)) {
                        diagnostics.recordClear(CanonicalTransitionCacheInvalidationReason.BRANCH_LENGTH_CHANGED);
                    }
                    if (entry.transition == null) {
                        entry.transition = new CanonicalGaussianTransition(getDimension());
                    }
                    ensureMomentBuffers(childNodeIndex);
                    transitionBuilder.fillCanonicalTransitionAndMoments(
                            snapshot,
                            effectiveBranchLength,
                            entry.transition,
                            transitionMomentMatrices[childNodeIndex],
                            transitionMomentCovariances[childNodeIndex]);
                    entry.processIndex = 0;
                    entry.updateSnapshot(childNodeIndex, effectiveBranchLength, entry.transition, null);
                    entry.markTransitionCurrent(effectiveBranchLength, epoch);
                } else {
                    diagnostics.recordHit();
                }
            }
        } else {
            diagnostics.recordHit();
        }
        return entry;
    }

    private void ensureMomentBuffers(final int childNodeIndex) {
        final int d2 = getDimension() * getDimension();
        if (transitionMomentMatrices[childNodeIndex] == null) {
            transitionMomentMatrices[childNodeIndex] = new double[d2];
        }
        if (transitionMomentCovariances[childNodeIndex] == null) {
            transitionMomentCovariances[childNodeIndex] = new double[d2];
        }
    }

    private void markDirty(final CanonicalTransitionCacheInvalidationReason reason) {
        dirty = true;
        if (reason == CanonicalTransitionCacheInvalidationReason.MODEL_CHANGED
                && dirtyReason != CanonicalTransitionCacheInvalidationReason.MODEL_CHANGED) {
            return;
        }
        dirtyReason = reason;
    }

    private void clearEntries(final CanonicalTransitionCacheInvalidationReason reason) {
        diagnostics.recordClear(reason);
        epoch.advance(reason);
        for (CanonicalTransitionCacheEntry entry : entries) {
            if (entry.invalidate(reason)) {
                diagnostics.recordPreparedCovarianceInvalidation();
            }
        }
    }

    private static void setIdentity(final MatrixParameter matrix) {
        for (int i = 0; i < matrix.getRowDimension(); ++i) {
            for (int j = 0; j < matrix.getColumnDimension(); ++j) {
                matrix.setParameterValueQuietly(i, j, i == j ? 1.0 : 0.0);
            }
        }
        matrix.fireParameterChangedEvent();
    }

    private void copyBlockBasisDiffusionCovariance() {
        final double[] qTilde = snapshot.getQDBasis().data;
        final int d = diffusionModel.getDimension();
        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d; ++j) {
                blockBasisDiffusionCovariance.setParameterValueQuietly(i, j, qTilde[i * d + j]);
            }
        }
    }
}

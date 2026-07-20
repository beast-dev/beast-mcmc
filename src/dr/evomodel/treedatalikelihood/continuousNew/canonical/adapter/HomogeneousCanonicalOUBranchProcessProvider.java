package dr.evomodel.treedatalikelihood.continuousNew.canonical.adapter;

import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.continuous.MultivariateDiffusionModel;
import dr.evomodel.continuous.MultivariateElasticModel;
import dr.evomodel.continuous.ou.OUProcessModel;
import dr.evomodel.continuous.ou.SelectionMatrixParameterization;
import dr.evomodel.treedatalikelihood.continuousNew.ContinuousBranchTimeModel;
import dr.evomodel.treedatalikelihood.continuous.ContinuousRateTransformation;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.CanonicalNativeBranchMessageProvider;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.message.CanonicalGaussianTransition;
import dr.inference.model.AbstractModel;
import dr.inference.model.MatrixParameter;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

/**
 * One OU process model shared by every non-root tree branch.
 */
final class HomogeneousCanonicalOUBranchProcessProvider implements CanonicalOUBranchProcessProvider {

    private final Tree tree;
    private final MultivariateDiffusionModel diffusionModel;
    private final MultivariateElasticModel elasticModel;
    private final MatrixParameterInterface driftMatrix;
    private final Parameter stationaryMean;
    private final MatrixParameter diffusionCovariance;
    private final CanonicalOUDiffusionSnapshot diffusionSnapshot;
    private final CanonicalBranchLengthSnapshot branchLengths;
    private final OUProcessModel processModel;

    HomogeneousCanonicalOUBranchProcessProvider(final Tree tree,
                                                final MultivariateElasticModel elasticModel,
                                                final MultivariateDiffusionModel diffusionModel,
                                                final Parameter stationaryMean,
                                                final BranchRateModel rateModel,
                                                final ContinuousRateTransformation rateTransformation) {
        this.tree = tree;
        this.elasticModel = elasticModel;
        this.diffusionModel = diffusionModel;
        this.driftMatrix = elasticModel.getStrengthOfSelectionMatrixParameter();
        this.stationaryMean = stationaryMean;
        final int dimension = driftMatrix.getRowDimension();
        this.diffusionCovariance = new MatrixParameter("canonicalOuProvider.diffusion", dimension, dimension);
        this.diffusionSnapshot = new CanonicalOUDiffusionSnapshot(
                dimension,
                diffusionModel,
                diffusionCovariance);
        this.branchLengths = new CanonicalBranchLengthSnapshot(
                tree,
                rateModel,
                rateTransformation);

        final MatrixParameter initialCovariance = new MatrixParameter("canonicalOuProvider.initial", dimension, dimension);
        setIdentity(initialCovariance);

        this.processModel = new OUProcessModel(
                "canonicalOuProvider.process",
                dimension,
                driftMatrix,
                diffusionCovariance,
                stationaryMean,
                initialCovariance);
    }

    @Override
    public int getDimension() {
        return driftMatrix.getRowDimension();
    }

    @Override
    public int getNodeCount() {
        return tree.getNodeCount();
    }

    @Override
    public int getProcessModelCount() {
        return 1;
    }

    @Override
    public OUProcessModel getProcessModel(final int processIndex) {
        if (processIndex != 0) {
            throw new IllegalArgumentException("Homogeneous OU process index must be 0: " + processIndex);
        }
        return processModel;
    }

    @Override
    public int getSegmentCount(final int childNodeIndex) {
        return branchLengths.getSegmentCount(childNodeIndex);
    }

    @Override
    public int getSegmentProcessIndex(final int childNodeIndex, final int segmentIndex) {
        branchLengths.getSegmentProcessIndex(childNodeIndex, segmentIndex);
        return 0;
    }

    @Override
    public double getSegmentEffectiveTime(final int childNodeIndex, final int segmentIndex) {
        return branchLengths.getSegmentEffectiveTime(childNodeIndex, segmentIndex);
    }

    @Override
    public OUProcessModel getProcessModel() {
        return processModel;
    }

    @Override
    public ContinuousBranchTimeModel getBranchTimeModel() {
        return branchLengths;
    }

    @Override
    public double getEffectiveBranchLength(final int childNodeIndex) {
        return branchLengths.getEffectiveBranchLength(childNodeIndex);
    }

    @Override
    public void fillTraitCovarianceFlat(final double[] out) {
        diffusionSnapshot.fillCovarianceFlat(out);
    }

    @Override
    public void fillCanonicalTransitionForLength(final double branchLength,
                                                 final CanonicalGaussianTransition out) {
        processModel.fillCanonicalTransition(branchLength, out);
    }

    @Override
    public void fillTransitionMatrixFlat(final double branchLength, final double[] out) {
        processModel.fillTransitionMatrixFlat(branchLength, out);
    }

    @Override
    public void fillTransitionOffset(final double branchLength, final double[] out) {
        processModel.fillTransitionOffset(branchLength, out);
    }

    @Override
    public void fillTransitionCovarianceFlat(final double branchLength, final double[] out) {
        processModel.fillTransitionCovarianceFlat(branchLength, out);
    }

    @Override
    public void fillTransitionMomentsFlat(final double branchLength,
                                          final double[] transitionMatrixOut,
                                          final double[] transitionOffsetOut,
                                          final double[] transitionCovarianceOut) {
        processModel.fillTransitionMomentsFlat(
                branchLength,
                transitionMatrixOut,
                transitionOffsetOut,
                transitionCovarianceOut);
    }

    @Override
    public CanonicalNativeBranchMessageProvider createNativeBranchMessageProvider(final int nodeCount) {
        final SelectionMatrixParameterization selection =
                processModel.getSelectionMatrixParameterization();
        return CanonicalNativeBranchMessageAdapterFactory.create(
                selection,
                nodeCount,
                diffusionCovariance,
                stationaryMean,
                branchLengths);
    }

    @Override
    public void registerWith(final AbstractModel owner) {
        owner.addModel(diffusionModel);
        owner.addModel(processModel);
        owner.addVariable(driftMatrix);
        owner.addVariable(stationaryMean);
    }

    @Override
    public void refreshSnapshot() {
        diffusionSnapshot.refresh();
        processModel.fireModelChanged();
    }

    @Override
    public void beginPhase() {
        branchLengths.beginPhase();
    }

    @Override
    public void endPhase() {
        branchLengths.endPhase();
    }

    @Override
    public CanonicalTransitionCacheInvalidationReason classifyModelChange(final Model model) {
        if (model == diffusionModel) {
            return CanonicalTransitionCacheInvalidationReason.DIFFUSION_CHANGED;
        }
        if (model == elasticModel) {
            return CanonicalTransitionCacheInvalidationReason.SELECTION_CHANGED;
        }
        return CanonicalTransitionCacheInvalidationReason.MODEL_CHANGED;
    }

    @Override
    public CanonicalTransitionCacheInvalidationReason classifyVariableChange(final Variable variable,
                                                                            final Parameter.ChangeType type) {
        if (variable == driftMatrix) {
            return CanonicalTransitionCacheInvalidationReason.SELECTION_CHANGED;
        }
        if (variable == stationaryMean) {
            return CanonicalTransitionCacheInvalidationReason.STATIONARY_MEAN_CHANGED;
        }
        return CanonicalTransitionCacheInvalidationReason.MODEL_CHANGED;
    }

    private static void setIdentity(final MatrixParameter matrix) {
        for (int i = 0; i < matrix.getRowDimension(); i++) {
            for (int j = 0; j < matrix.getColumnDimension(); j++) {
                matrix.setParameterValueQuietly(i, j, i == j ? 1.0 : 0.0);
            }
        }
        matrix.fireParameterChangedEvent();
    }
}

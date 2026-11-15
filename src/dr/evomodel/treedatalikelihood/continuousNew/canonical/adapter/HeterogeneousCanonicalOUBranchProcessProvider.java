package dr.evomodel.treedatalikelihood.continuousNew.canonical.adapter;

import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.continuous.MultivariateDiffusionModel;
import dr.evomodel.continuous.MultivariateElasticModel;
import dr.evomodel.continuous.ou.OUProcessModel;
import dr.evomodel.continuous.ou.SelectionMatrixParameterization;
import dr.evomodel.treedatalikelihood.continuousNew.CanonicalOUBranchProcessSpecification;
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
 * One-segment branch-process provider with branch-specific OU selection models.
 *
 * <p>The first heterogeneous canonical OU path keeps diffusion and stationary mean
 * shared across processes. Branches select a process by child node index through a
 * mutable integer-valued parameter, giving likelihood and gradients a common source
 * of branch-local process identity.</p>
 */
final class HeterogeneousCanonicalOUBranchProcessProvider implements CanonicalOUBranchProcessProvider {

    private static final double INTEGER_TOLERANCE = 1.0e-9;

    private final Tree tree;
    private final MultivariateDiffusionModel diffusionModel;
    private final MultivariateElasticModel[] elasticModels;
    private final MatrixParameterInterface[] driftMatrices;
    private final Parameter stationaryMean;
    private final Parameter branchProcessIndexParameter;
    private final int[] assignmentIndexByNode;
    private final MatrixParameter diffusionCovariance;
    private final CanonicalOUDiffusionSnapshot diffusionSnapshot;
    private final CanonicalBranchLengthSnapshot branchLengths;
    private final OUProcessModel[] processModels;
    private final int dimension;

    HeterogeneousCanonicalOUBranchProcessProvider(
            final Tree tree,
            final CanonicalOUBranchProcessSpecification specification,
            final MultivariateDiffusionModel diffusionModel,
            final Parameter stationaryMean,
            final BranchRateModel rateModel,
            final ContinuousRateTransformation rateTransformation) {
        this.tree = tree;
        this.diffusionModel = diffusionModel;
        this.stationaryMean = stationaryMean;
        this.branchProcessIndexParameter = specification.getBranchProcessIndexParameter();

        final int processCount = specification.getProcessCount();
        this.elasticModels = new MultivariateElasticModel[processCount];
        this.driftMatrices = new MatrixParameterInterface[processCount];
        for (int i = 0; i < processCount; ++i) {
            this.elasticModels[i] = specification.getElasticModel(i);
            this.driftMatrices[i] = elasticModels[i].getStrengthOfSelectionMatrixParameter();
        }

        this.dimension = requireConsistentSelectionDimension(driftMatrices);
        this.diffusionCovariance =
                new MatrixParameter("canonicalOuHeterogeneousProvider.diffusion", dimension, dimension);
        this.diffusionSnapshot = new CanonicalOUDiffusionSnapshot(
                dimension,
                diffusionModel,
                diffusionCovariance);
        this.branchLengths = new CanonicalBranchLengthSnapshot(
                tree,
                rateModel,
                rateTransformation);
        this.assignmentIndexByNode = buildAssignmentIndexByNode(tree, branchProcessIndexParameter);
        this.processModels = new OUProcessModel[processCount];

        for (int i = 0; i < processCount; ++i) {
            final MatrixParameter initialCovariance =
                    new MatrixParameter("canonicalOuHeterogeneousProvider.initial." + i, dimension, dimension);
            setIdentity(initialCovariance);
            processModels[i] = new OUProcessModel(
                    "canonicalOuHeterogeneousProvider.process." + i,
                    dimension,
                    driftMatrices[i],
                    diffusionCovariance,
                    stationaryMean,
                    initialCovariance);
        }
    }

    @Override
    public int getDimension() {
        return dimension;
    }

    @Override
    public int getNodeCount() {
        return tree.getNodeCount();
    }

    @Override
    public int getProcessModelCount() {
        return processModels.length;
    }

    @Override
    public OUProcessModel getProcessModel(final int processIndex) {
        return processModels[processIndex];
    }

    @Override
    public int getSegmentCount(final int childNodeIndex) {
        return branchLengths.getSegmentCount(childNodeIndex);
    }

    @Override
    public int getSegmentProcessIndex(final int childNodeIndex, final int segmentIndex) {
        branchLengths.getSegmentProcessIndex(childNodeIndex, segmentIndex);
        return processIndexForChild(childNodeIndex);
    }

    @Override
    public double getSegmentEffectiveTime(final int childNodeIndex, final int segmentIndex) {
        return branchLengths.getSegmentEffectiveTime(childNodeIndex, segmentIndex);
    }

    @Override
    public OUProcessModel getProcessModel() {
        return processModels[0];
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
        processModels[0].fillCanonicalTransition(branchLength, out);
    }

    @Override
    public void fillCanonicalTransitionForBranchLength(final int childNodeIndex,
                                                       final double branchLength,
                                                       final CanonicalGaussianTransition out) {
        processModels[processIndexForChild(childNodeIndex)].fillCanonicalTransition(branchLength, out);
    }

    @Override
    public void fillTransitionMatrixFlat(final double branchLength, final double[] out) {
        processModels[0].fillTransitionMatrixFlat(branchLength, out);
    }

    @Override
    public void fillTransitionMatrixFlatForBranchLength(final int childNodeIndex,
                                                        final double branchLength,
                                                        final double[] out) {
        processModels[processIndexForChild(childNodeIndex)].fillTransitionMatrixFlat(branchLength, out);
    }

    @Override
    public void fillTransitionOffset(final double branchLength, final double[] out) {
        processModels[0].fillTransitionOffset(branchLength, out);
    }

    @Override
    public void fillTransitionOffsetForBranchLength(final int childNodeIndex,
                                                    final double branchLength,
                                                    final double[] out) {
        processModels[processIndexForChild(childNodeIndex)].fillTransitionOffset(branchLength, out);
    }

    @Override
    public void fillTransitionCovarianceFlat(final double branchLength, final double[] out) {
        processModels[0].fillTransitionCovarianceFlat(branchLength, out);
    }

    @Override
    public void fillTransitionCovarianceFlatForBranchLength(final int childNodeIndex,
                                                            final double branchLength,
                                                            final double[] out) {
        processModels[processIndexForChild(childNodeIndex)].fillTransitionCovarianceFlat(branchLength, out);
    }

    @Override
    public void fillTransitionMomentsFlat(final double branchLength,
                                          final double[] transitionMatrixOut,
                                          final double[] transitionOffsetOut,
                                          final double[] transitionCovarianceOut) {
        processModels[0].fillTransitionMomentsFlat(
                branchLength,
                transitionMatrixOut,
                transitionOffsetOut,
                transitionCovarianceOut);
    }

    @Override
    public void fillTransitionMomentsFlatForBranchLength(final int childNodeIndex,
                                                         final double branchLength,
                                                         final double[] transitionMatrixOut,
                                                         final double[] transitionOffsetOut,
                                                         final double[] transitionCovarianceOut) {
        processModels[processIndexForChild(childNodeIndex)].fillTransitionMomentsFlat(
                branchLength,
                transitionMatrixOut,
                transitionOffsetOut,
                transitionCovarianceOut);
    }

    @Override
    public CanonicalNativeBranchMessageProvider createNativeBranchMessageProvider(final int nodeCount) {
        if (processModels.length != 1) {
            return null;
        }
        final SelectionMatrixParameterization selection =
                processModels[0].getSelectionMatrixParameterization();
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
        owner.addVariable(stationaryMean);
        owner.addVariable(branchProcessIndexParameter);
        for (int i = 0; i < processModels.length; ++i) {
            owner.addModel(elasticModels[i]);
            owner.addModel(processModels[i]);
            owner.addVariable(driftMatrices[i]);
        }
    }

    @Override
    public void refreshSnapshot() {
        diffusionSnapshot.refresh();
        for (int i = 0; i < processModels.length; ++i) {
            processModels[i].fireModelChanged();
        }
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
        for (int i = 0; i < processModels.length; ++i) {
            if (model == elasticModels[i]) {
                return CanonicalTransitionCacheInvalidationReason.SELECTION_CHANGED;
            }
            if (model == processModels[i]) {
                return CanonicalTransitionCacheInvalidationReason.MODEL_CHANGED;
            }
        }
        return CanonicalTransitionCacheInvalidationReason.MODEL_CHANGED;
    }

    @Override
    public CanonicalTransitionCacheInvalidationReason classifyVariableChange(final Variable variable,
                                                                            final Parameter.ChangeType type) {
        if (variable == stationaryMean) {
            return CanonicalTransitionCacheInvalidationReason.STATIONARY_MEAN_CHANGED;
        }
        if (variable == branchProcessIndexParameter) {
            return CanonicalTransitionCacheInvalidationReason.MODEL_CHANGED;
        }
        for (int i = 0; i < driftMatrices.length; ++i) {
            if (variable == driftMatrices[i]) {
                return CanonicalTransitionCacheInvalidationReason.SELECTION_CHANGED;
            }
        }
        return CanonicalTransitionCacheInvalidationReason.MODEL_CHANGED;
    }

    private int processIndexForChild(final int childNodeIndex) {
        final int parameterIndex = assignmentIndexByNode[childNodeIndex];
        if (parameterIndex < 0) {
            return 0;
        }
        final double raw = branchProcessIndexParameter.getParameterValue(parameterIndex);
        final int processIndex = (int) Math.rint(raw);
        if (Math.abs(raw - processIndex) > INTEGER_TOLERANCE) {
            throw new IllegalArgumentException(
                    "OU branch-process index for child node " + childNodeIndex
                            + " must be integer-valued, found " + raw + ".");
        }
        if (processIndex < 0 || processIndex >= processModels.length) {
            throw new IndexOutOfBoundsException(
                    "OU branch-process index for child node " + childNodeIndex
                            + " is " + processIndex + ", but valid range is [0,"
                            + (processModels.length - 1) + "].");
        }
        return processIndex;
    }

    private static int requireConsistentSelectionDimension(
            final MatrixParameterInterface[] driftMatrices) {
        final int dimension = driftMatrices[0].getRowDimension();
        if (dimension < 1 || driftMatrices[0].getColumnDimension() != dimension) {
            throw new IllegalArgumentException("OU selection matrix must be square.");
        }
        for (int i = 1; i < driftMatrices.length; ++i) {
            if (driftMatrices[i].getRowDimension() != dimension
                    || driftMatrices[i].getColumnDimension() != dimension) {
                throw new IllegalArgumentException(
                        "All OU branch-process selection matrices must have dimension "
                                + dimension + " x " + dimension + ".");
            }
        }
        return dimension;
    }

    private static int[] buildAssignmentIndexByNode(final Tree tree,
                                                    final Parameter branchProcessIndexParameter) {
        final int nodeCount = tree.getNodeCount();
        final int dimension = branchProcessIndexParameter.getDimension();
        final int rootIndex = tree.getRoot().getNumber();
        final int[] indexByNode = new int[nodeCount];

        if (dimension == nodeCount) {
            for (int i = 0; i < nodeCount; ++i) {
                indexByNode[i] = i;
            }
            indexByNode[rootIndex] = -1;
            return indexByNode;
        }

        if (dimension == nodeCount - 1) {
            int parameterIndex = 0;
            for (int nodeIndex = 0; nodeIndex < nodeCount; ++nodeIndex) {
                if (nodeIndex == rootIndex) {
                    indexByNode[nodeIndex] = -1;
                } else {
                    indexByNode[nodeIndex] = parameterIndex++;
                }
            }
            return indexByNode;
        }

        throw new IllegalArgumentException(
                "OU branch-process index parameter dimension must be nodeCount ("
                        + nodeCount + ") or nodeCount - 1 (" + (nodeCount - 1)
                        + "), found " + dimension + ".");
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

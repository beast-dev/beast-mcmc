package dr.inference.timeseries.gaussian;

import dr.evomodel.continuous.ou.OUProcessModel;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianBranchTransitionKernel;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.GaussianBranchTransitionKernel;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.inference.timeseries.core.LatentProcessModel;
import dr.inference.timeseries.representation.GaussianComputationMode;
import dr.inference.timeseries.representation.CachedGaussianTransitionRepresentation;
import dr.inference.timeseries.representation.GaussianTransitionRepresentation;
import dr.inference.timeseries.representation.KernelBackedGaussianTransitionRepresentation;
import dr.inference.timeseries.representation.RepresentableProcess;

import java.util.IdentityHashMap;

/**
 * Time-series representation adapter for the exact OU process model.
 *
 * <p>The OU process itself owns the branch-length transition mathematics. This
 * adapter owns the time-grid representation contract used by Kalman engines.</p>
 */
public final class OUTimeSeriesProcessAdapter extends AbstractModel
        implements LatentProcessModel, RepresentableProcess {

    private static final IdentityHashMap<OUProcessModel, KernelBackedGaussianTransitionRepresentation>
            SHARED_TRANSITION_REPRESENTATIONS =
            new IdentityHashMap<OUProcessModel, KernelBackedGaussianTransitionRepresentation>();

    private final OUProcessModel processModel;
    private final GaussianComputationMode defaultComputationMode;
    private final GaussianTransitionRepresentation transitionRepresentation;

    public OUTimeSeriesProcessAdapter(final OUProcessModel processModel) {
        this(processModel, GaussianComputationMode.EXPECTATION);
    }

    public OUTimeSeriesProcessAdapter(final OUProcessModel processModel,
                                      final GaussianComputationMode defaultComputationMode) {
        super(processModel == null ? "ouTimeSeriesProcessAdapter" : processModel.getId() + ".timeSeries");
        if (processModel == null) {
            throw new IllegalArgumentException("processModel must not be null");
        }
        if (defaultComputationMode == null) {
            throw new IllegalArgumentException("defaultComputationMode must not be null");
        }
        this.processModel = processModel;
        this.defaultComputationMode = defaultComputationMode;
        this.transitionRepresentation = sharedTransitionRepresentationFor(processModel);
        addModel(processModel);
    }

    private static synchronized KernelBackedGaussianTransitionRepresentation sharedTransitionRepresentationFor(
            final OUProcessModel processModel) {
        KernelBackedGaussianTransitionRepresentation representation =
                SHARED_TRANSITION_REPRESENTATIONS.get(processModel);
        if (representation == null) {
            representation = new KernelBackedGaussianTransitionRepresentation(processModel);
            SHARED_TRANSITION_REPRESENTATIONS.put(processModel, representation);
        }
        return representation;
    }

    public OUProcessModel getProcessModel() {
        return processModel;
    }

    public GaussianComputationMode getDefaultComputationMode() {
        return defaultComputationMode;
    }

    @Override
    public int getStateDimension() {
        return processModel.getStateDimension();
    }

    @Override
    public <T> boolean supportsRepresentation(final Class<T> representationClass) {
        return representationClass.isAssignableFrom(GaussianTransitionRepresentation.class)
                || representationClass.isAssignableFrom(GaussianBranchTransitionKernel.class)
                || representationClass.isAssignableFrom(CanonicalGaussianBranchTransitionKernel.class);
    }

    @Override
    @SuppressWarnings("unchecked")
    public <T> T getRepresentation(final Class<T> representationClass) {
        if (!supportsRepresentation(representationClass)) {
            throw new IllegalArgumentException("Unsupported representation: " + representationClass.getName());
        }
        if (representationClass.isAssignableFrom(GaussianTransitionRepresentation.class)) {
            return (T) transitionRepresentation;
        }
        return (T) processModel;
    }

    @Override
    protected void handleModelChangedEvent(final Model model, final Object object, final int index) {
        if (transitionRepresentation instanceof CachedGaussianTransitionRepresentation) {
            ((CachedGaussianTransitionRepresentation) transitionRepresentation).makeDirty();
        }
        fireModelChanged();
    }

    @Override
    protected void handleVariableChangedEvent(final Variable variable,
                                              final int index,
                                              final Parameter.ChangeType type) {
        if (transitionRepresentation instanceof CachedGaussianTransitionRepresentation) {
            ((CachedGaussianTransitionRepresentation) transitionRepresentation).makeDirty();
        }
        fireModelChanged();
    }

    @Override
    protected void storeState() {
        // no-op
    }

    @Override
    protected void restoreState() {
        if (transitionRepresentation instanceof CachedGaussianTransitionRepresentation) {
            ((CachedGaussianTransitionRepresentation) transitionRepresentation).makeDirty();
        }
    }

    @Override
    protected void acceptState() {
        // no-op
    }
}

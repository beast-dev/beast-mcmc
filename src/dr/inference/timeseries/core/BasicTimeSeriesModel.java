package dr.inference.timeseries.core;

import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

/**
 * Minimal concrete implementation of {@link TimeSeriesModel}.
 */
public class BasicTimeSeriesModel extends AbstractModel implements TimeSeriesModel {

    private final LatentProcessModel latentProcessModel;
    private final ObservationModel observationModel;
    private final TimeGrid timeGrid;

    public BasicTimeSeriesModel(final String name,
                                final LatentProcessModel latentProcessModel,
                                final ObservationModel observationModel,
                                final TimeGrid timeGrid) {
        super(name);
        if (latentProcessModel == null) {
            throw new IllegalArgumentException("latentProcessModel must not be null");
        }
        if (observationModel == null) {
            throw new IllegalArgumentException("observationModel must not be null");
        }
        if (timeGrid == null) {
            throw new IllegalArgumentException("timeGrid must not be null");
        }
        this.latentProcessModel = latentProcessModel;
        this.observationModel = observationModel;
        this.timeGrid = timeGrid;

        addModel(latentProcessModel);
        addModel(observationModel);
    }

    @Override
    public LatentProcessModel getLatentProcessModel() {
        return latentProcessModel;
    }

    @Override
    public ObservationModel getObservationModel() {
        return observationModel;
    }

    @Override
    public TimeGrid getTimeGrid() {
        return timeGrid;
    }

    @Override
    protected void handleModelChangedEvent(final Model model, final Object object, final int index) {
        fireModelChanged();
    }

    @Override
    protected void handleVariableChangedEvent(final Variable variable,
                                              final int index,
                                              final Parameter.ChangeType type) {
        fireModelChanged();
    }

    @Override
    protected void storeState() {
        // no-op
    }

    @Override
    protected void restoreState() {
        // no-op
    }

    @Override
    protected void acceptState() {
        // no-op
    }
}

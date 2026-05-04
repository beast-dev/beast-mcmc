package dr.inference.timeseries.core;

import dr.inference.model.Model;

/**
 * Semantic container tying together a latent process, an observation model,
 * and a time grid.
 */
public interface TimeSeriesModel extends Model {

    LatentProcessModel getLatentProcessModel();

    ObservationModel getObservationModel();

    TimeGrid getTimeGrid();
}

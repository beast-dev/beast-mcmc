package dr.inference.timeseries.core;

import dr.inference.model.Model;

/**
 * Semantic interface for an observation model attached to a latent process.
 */
public interface ObservationModel extends Model {

    /**
     * @return dimension of the observed vector at each observation time.
     */
    int getObservationDimension();
}

package dr.inference.timeseries.core;

import dr.inference.model.Model;

/**
 * Semantic interface for a latent time-series process.
 *
 * This interface intentionally stays small: it describes what the latent process is,
 * not how it is evaluated numerically.
 */
public interface LatentProcessModel extends Model {

    /**
     * @return dimension of the latent state at each time point.
     */
    int getStateDimension();
}

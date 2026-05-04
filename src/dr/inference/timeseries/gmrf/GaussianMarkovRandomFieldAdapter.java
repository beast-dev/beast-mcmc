package dr.inference.timeseries.gmrf;

import dr.inference.timeseries.representation.SparsePrecisionRepresentation;
import dr.math.distributions.GaussianMarkovRandomField;

/**
 * Adapter exposing the current BEAST GaussianMarkovRandomField implementation as a
 * sparse precision representation in the new time-series architecture.
 */
public class GaussianMarkovRandomFieldAdapter implements SparsePrecisionRepresentation {

    private final GaussianMarkovRandomField gmrf;

    public GaussianMarkovRandomFieldAdapter(final GaussianMarkovRandomField gmrf) {
        if (gmrf == null) {
            throw new IllegalArgumentException("gmrf must not be null");
        }
        this.gmrf = gmrf;
    }

    @Override
    public int getDimension() {
        return gmrf.getDimension();
    }

    @Override
    public int getBandwidth() {
        return 1;
    }

    @Override
    public double getEntry(final int row, final int col) {
        throw new UnsupportedOperationException(
                "The wrapped GaussianMarkovRandomField does not expose raw precision entries through this adapter; " +
                "use the native GMRF API for entry-level access.");
    }

    @Override
    public double getLogDeterminant() {
        return gmrf.getLogDeterminant();
    }

    @Override
    public void updateStructure() {
        // The wrapped GMRF already manages its own internal cache invalidation.
    }

    public GaussianMarkovRandomField getWrappedGMRF() {
        return gmrf;
    }
}

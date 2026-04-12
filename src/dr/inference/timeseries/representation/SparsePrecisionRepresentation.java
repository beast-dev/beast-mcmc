package dr.inference.timeseries.representation;

/**
 * Computational representation for models naturally expressed by sparse precision matrices.
 */
public interface SparsePrecisionRepresentation {

    int getDimension();

    int getBandwidth();

    double getEntry(int row, int col);

    double getLogDeterminant();

    void updateStructure();
}

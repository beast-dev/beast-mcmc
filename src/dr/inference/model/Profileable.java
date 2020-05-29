package dr.inference.model;

/**
 * An interface to add into a likelihood to allow for granular profiling of performance.
 * Currently, simply returns a count of the number of computations done (of whatever type is appropriate).
 * @author Andrew Rambaut
 * @version $Id$
 */
public interface Profileable {

    /**
     * returns the cumulative count of the number of individual calculations done by this likelihood.
     * @return
     */
    long getTotalCalculationCount();

}

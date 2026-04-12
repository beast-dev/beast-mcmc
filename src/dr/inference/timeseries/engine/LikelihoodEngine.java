package dr.inference.timeseries.engine;

/**
 * Internal likelihood engine.
 */
public interface LikelihoodEngine {

    double getLogLikelihood();

    void makeDirty();
}

package dr.inference.timeseries.engine.gmrf;

import dr.inference.model.Parameter;
import dr.inference.timeseries.engine.LikelihoodEngine;
import dr.math.distributions.GaussianMarkovRandomField;

/**
 * Likelihood engine delegating to the existing GaussianMarkovRandomField implementation.
 */
public class SparsePrecisionLikelihoodEngine implements LikelihoodEngine {

    private final GaussianMarkovRandomField gmrf;
    private final Parameter fieldParameter;

    private boolean likelihoodKnown;
    private double logLikelihood;

    public SparsePrecisionLikelihoodEngine(final GaussianMarkovRandomField gmrf,
                                           final Parameter fieldParameter) {
        if (gmrf == null) {
            throw new IllegalArgumentException("gmrf must not be null");
        }
        if (fieldParameter == null) {
            throw new IllegalArgumentException("fieldParameter must not be null");
        }
        this.gmrf = gmrf;
        this.fieldParameter = fieldParameter;
    }

    @Override
    public double getLogLikelihood() {
        if (!likelihoodKnown) {
            logLikelihood = gmrf.logPdf(fieldParameter.getParameterValues());
            likelihoodKnown = true;
        }
        return logLikelihood;
    }

    @Override
    public void makeDirty() {
        likelihoodKnown = false;
    }
}

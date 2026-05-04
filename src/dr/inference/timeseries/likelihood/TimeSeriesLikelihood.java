package dr.inference.timeseries.likelihood;

import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.GradientProvider;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.inference.timeseries.core.TimeSeriesModel;
import dr.inference.timeseries.engine.GradientEngine;
import dr.inference.timeseries.engine.LikelihoodEngine;

/**
 * Top-level BEAST likelihood wrapper for time-series models.
 */
public class TimeSeriesLikelihood extends AbstractModelLikelihood implements TimeSeriesGradientSource {

    private final TimeSeriesModel model;
    private final LikelihoodEngine likelihoodEngine;
    private final GradientEngine gradientEngine;

    private boolean likelihoodKnown;
    private double logLikelihood;

    private boolean storedLikelihoodKnown;
    private double storedLogLikelihood;

    public TimeSeriesLikelihood(final String name,
                                final TimeSeriesModel model,
                                final LikelihoodEngine likelihoodEngine,
                                final GradientEngine gradientEngine) {
        super(name);
        if (model == null) {
            throw new IllegalArgumentException("model must not be null");
        }
        if (likelihoodEngine == null) {
            throw new IllegalArgumentException("likelihoodEngine must not be null");
        }
        if (gradientEngine == null) {
            throw new IllegalArgumentException("gradientEngine must not be null");
        }
        this.model = model;
        this.likelihoodEngine = likelihoodEngine;
        this.gradientEngine = gradientEngine;

        addModel(model);
        this.likelihoodKnown = false;
    }

    @Override
    public Model getModel() {
        return this;
    }

    public TimeSeriesModel getTimeSeriesModel() {
        return model;
    }

    @Override
    public double getLogLikelihood() {
        if (!likelihoodKnown) {
            logLikelihood = likelihoodEngine.getLogLikelihood();
            likelihoodKnown = true;
        }
        return logLikelihood;
    }

    @Override
    public void makeDirty() {
        likelihoodKnown = false;
        likelihoodEngine.makeDirty();
        gradientEngine.makeDirty();
    }

    @Override
    public GradientProvider getGradientWrt(final Parameter parameter) {
        if (!gradientEngine.supportsGradientWrt(parameter)) {
            throw new IllegalArgumentException("Unsupported gradient parameter: " + parameter.getParameterName());
        }

        return new GradientProvider() {
            @Override
            public int getDimension() {
                return parameter.getDimension();
            }

            @Override
            public double[] getGradientLogDensity(final Object x) {
                return gradientEngine.getGradientWrt(parameter);
            }
        };
    }

    @Override
    protected void handleModelChangedEvent(final Model model, final Object object, final int index) {
        makeDirty();
        fireModelChanged(object, index);
    }

    @Override
    protected void handleVariableChangedEvent(final Variable variable,
                                              final int index,
                                              final Parameter.ChangeType type) {
        makeDirty();
        fireModelChanged();
    }

    @Override
    protected void storeState() {
        storedLikelihoodKnown = likelihoodKnown;
        storedLogLikelihood = logLikelihood;
    }

    @Override
    protected void restoreState() {
        likelihoodKnown = storedLikelihoodKnown;
        logLikelihood = storedLogLikelihood;
        likelihoodEngine.makeDirty();
        gradientEngine.makeDirty();
    }

    @Override
    protected void acceptState() {
        // no-op
    }
}

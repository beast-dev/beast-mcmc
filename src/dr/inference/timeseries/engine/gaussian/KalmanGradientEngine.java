package dr.inference.timeseries.engine.gaussian;

import dr.inference.model.Parameter;
import dr.inference.timeseries.engine.GradientEngine;
import dr.inference.timeseries.gaussian.GaussianObservationModel;
import dr.evomodel.continuous.ou.OUProcessModel;

/**
 * Gradient engine for the Kalman likelihood path.
 *
 * This implementation provides a reliable central-difference fallback over BEAST parameters.
 * It is numerically expensive, but it lets the module run end-to-end immediately while keeping
 * the public architecture stable. It can later be replaced by an analytic backward pass.
 */
public class KalmanGradientEngine implements GradientEngine {

    private static final double DEFAULT_STEP_SCALE = 1E-6;
    private static final double MINIMUM_STEP = 1E-8;

    private final KalmanLikelihoodEngine likelihoodEngine;
    private final Parameter[] supportedParameters;

    public KalmanGradientEngine(final KalmanLikelihoodEngine likelihoodEngine,
                                final OUProcessModel processModel,
                                final GaussianObservationModel observationModel) {
        if (likelihoodEngine == null) {
            throw new IllegalArgumentException("likelihoodEngine must not be null");
        }
        if (processModel == null) {
            throw new IllegalArgumentException("processModel must not be null");
        }
        if (observationModel == null) {
            throw new IllegalArgumentException("observationModel must not be null");
        }
        this.likelihoodEngine = likelihoodEngine;
        this.supportedParameters = new Parameter[] {
                processModel.getDriftMatrix(),
                processModel.getDiffusionMatrix(),
                processModel.getStationaryMeanParameter(),
                processModel.getInitialCovarianceParameter(),
                observationModel.getDesignMatrix(),
                observationModel.getNoiseCovariance(),
                observationModel.getObservations()
        };
    }

    @Override
    public boolean supportsGradientWrt(final Parameter parameter) {
        if (parameter == null) {
            return false;
        }
        for (Parameter supportedParameter : supportedParameters) {
            if (parameter == supportedParameter) {
                return true;
            }
        }
        return false;
    }

    @Override
    public double[] getGradientWrt(final Parameter parameter) {
        if (!supportsGradientWrt(parameter)) {
            throw new IllegalArgumentException("Unsupported parameter for KalmanGradientEngine");
        }

        final int dimension = parameter.getDimension();
        final double[] gradient = new double[dimension];
        final double[] savedValues = parameter.getParameterValues();

        for (int i = 0; i < dimension; ++i) {
            final double originalValue = savedValues[i];
            final double step = stepSize(originalValue);

            try {
                setParameterValue(parameter, i, originalValue + step);
                final double plus = evaluateCurrentState();

                setParameterValue(parameter, i, originalValue - step);
                final double minus = evaluateCurrentState();

                gradient[i] = (plus - minus) / (2.0 * step);
            } finally {
                setParameterValue(parameter, i, originalValue);
            }
        }

        parameter.fireParameterChangedEvent();
        likelihoodEngine.makeDirty();
        return gradient;
    }

    @Override
    public void makeDirty() {
        likelihoodEngine.makeDirty();
    }

    private double evaluateCurrentState() {
        likelihoodEngine.makeDirty();
        return likelihoodEngine.getLogLikelihood();
    }

    private static double stepSize(final double value) {
        return Math.max(MINIMUM_STEP, DEFAULT_STEP_SCALE * Math.max(1.0, Math.abs(value)));
    }

    private static void setParameterValue(final Parameter parameter, final int index, final double value) {
        parameter.setParameterValueQuietly(index, value);
        parameter.fireParameterChangedEvent();
    }
}

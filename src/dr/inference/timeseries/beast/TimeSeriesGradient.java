package dr.inference.timeseries.beast;

import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.GradientProvider;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.inference.timeseries.likelihood.TimeSeriesLikelihood;
import dr.xml.Reportable;

/**
 * BEAST-facing HMC adapter for time-series likelihood gradients.
 */
public class TimeSeriesGradient implements GradientWrtParameterProvider, Reportable {

    private final TimeSeriesLikelihood likelihood;
    private final Parameter parameter;
    private final GradientProvider provider;

    public TimeSeriesGradient(final TimeSeriesLikelihood likelihood,
                              final Parameter parameter) {
        if (likelihood == null) {
            throw new IllegalArgumentException("likelihood must not be null");
        }
        if (parameter == null) {
            throw new IllegalArgumentException("parameter must not be null");
        }
        this.likelihood = likelihood;
        this.parameter = parameter;
        this.provider = likelihood.getGradientWrt(parameter);
    }

    @Override
    public String getReport() {
        return GradientWrtParameterProvider.getReportAndCheckForError(
                this,
                parameter.getBounds().getLowerLimit(0),
                parameter.getBounds().getUpperLimit(0),
                null
        );
    }

    @Override
    public Likelihood getLikelihood() {
        return likelihood;
    }

    @Override
    public Parameter getParameter() {
        return parameter;
    }

    @Override
    public int getDimension() {
        return parameter.getDimension();
    }

    @Override
    public double[] getGradientLogDensity() {
        return provider.getGradientLogDensity(parameter.getParameterValues());
    }
}

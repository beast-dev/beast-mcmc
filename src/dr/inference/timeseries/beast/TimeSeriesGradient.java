package dr.inference.timeseries.beast;

import dr.inference.hmc.BatchGradientWrtParameterProvider;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.GradientProvider;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.inference.timeseries.likelihood.ParallelTimeSeriesLikelihood;
import dr.inference.timeseries.likelihood.TimeSeriesGradientSource;
import dr.xml.Reportable;

import java.util.ArrayList;
import java.util.List;

/**
 * BEAST-facing HMC adapter for time-series likelihood gradients.
 */
public class TimeSeriesGradient implements BatchGradientWrtParameterProvider, Reportable {

    private final TimeSeriesGradientSource gradientSource;
    private final Likelihood likelihood;
    private final Parameter parameter;
    private final GradientProvider provider;

    public TimeSeriesGradient(final TimeSeriesGradientSource likelihood,
                              final Parameter parameter) {
        if (likelihood == null) {
            throw new IllegalArgumentException("likelihood must not be null");
        }
        if (!(likelihood instanceof Likelihood)) {
            throw new IllegalArgumentException("likelihood must also implement BEAST Likelihood");
        }
        if (parameter == null) {
            throw new IllegalArgumentException("parameter must not be null");
        }
        this.gradientSource = likelihood;
        this.likelihood = (Likelihood) likelihood;
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
        return provider.getGradientLogDensity(null);
    }

    @Override
    public Object getBatchGradientKey() {
        return gradientSource instanceof ParallelTimeSeriesLikelihood ? gradientSource : null;
    }

    @Override
    public double[][] getGradientLogDensityBatch(final List<BatchGradientWrtParameterProvider> providers) {
        if (!(gradientSource instanceof ParallelTimeSeriesLikelihood)) {
            final double[][] gradients = new double[providers.size()][];
            for (int i = 0; i < providers.size(); ++i) {
                gradients[i] = providers.get(i).getGradientLogDensity();
            }
            return gradients;
        }
        final List<Parameter> parameters = new ArrayList<Parameter>(providers.size());
        for (final BatchGradientWrtParameterProvider provider : providers) {
            if (!(provider instanceof TimeSeriesGradient)
                    || ((TimeSeriesGradient) provider).gradientSource != gradientSource) {
                throw new IllegalArgumentException("All batched time-series gradients must share one likelihood.");
            }
            parameters.add(provider.getParameter());
        }
        return ((ParallelTimeSeriesLikelihood) gradientSource).computeGradients(parameters);
    }
}

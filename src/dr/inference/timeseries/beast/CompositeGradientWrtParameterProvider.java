package dr.inference.timeseries.beast;

import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Concatenates several gradient providers that share one likelihood.
 */
public final class CompositeGradientWrtParameterProvider implements GradientWrtParameterProvider {

    private final List<GradientWrtParameterProvider> providers;
    private final Likelihood likelihood;
    private final CompoundParameter parameter;
    private final int dimension;

    public CompositeGradientWrtParameterProvider(final String name,
                                                 final List<GradientWrtParameterProvider> providers) {
        if (providers == null || providers.isEmpty()) {
            throw new IllegalArgumentException("providers must not be empty");
        }
        this.providers = Collections.unmodifiableList(
                new ArrayList<GradientWrtParameterProvider>(providers));
        this.likelihood = providers.get(0).getLikelihood();
        this.parameter = new CompoundParameter(name == null ? "compositeGradientParameter" : name);

        int totalDimension = 0;
        for (final GradientWrtParameterProvider provider : providers) {
            if (provider == null) {
                throw new IllegalArgumentException("providers must not contain null");
            }
            if (provider.getLikelihood() != likelihood) {
                throw new IllegalArgumentException("all providers must share one likelihood");
            }
            parameter.addParameter(provider.getParameter());
            totalDimension += provider.getDimension();
        }
        this.dimension = totalDimension;
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
        return dimension;
    }

    @Override
    public double[] getGradientLogDensity() {
        final double[] gradient = new double[dimension];
        int offset = 0;
        for (final GradientWrtParameterProvider provider : providers) {
            final double[] part = provider.getGradientLogDensity();
            if (part.length != provider.getDimension()) {
                throw new IllegalStateException(
                        "gradient length mismatch for " + provider.getParameter().getParameterName());
            }
            System.arraycopy(part, 0, gradient, offset, part.length);
            offset += part.length;
        }
        return gradient;
    }
}

package dr.inference.hmc;

import java.util.List;

/**
 * Optional fast path for providers that can compute several parameter gradients
 * together while sharing upstream work.
 */
public interface BatchGradientWrtParameterProvider extends GradientWrtParameterProvider {

    Object getBatchGradientKey();

    double[][] getGradientLogDensityBatch(List<BatchGradientWrtParameterProvider> providers);
}

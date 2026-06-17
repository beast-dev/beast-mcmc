package dr.math.distributions.gp;

import dr.inference.distribution.RandomField;
import dr.inference.distribution.RandomFieldGradient;
import dr.inference.model.DesignMatrix;
import dr.inference.model.GradientProvider;
import dr.inference.model.Parameter;

/**
 * Gradient provider for Gaussian process random fields.
 */
public class GaussianProcessFieldGradient extends RandomFieldGradient {

    public GaussianProcessFieldGradient(RandomField randomField,
                                        Parameter parameter) {
        super(randomField, parameter);
        getGaussianProcessField(randomField);
    }

    static GradientProvider createProvider(final GaussianProcessField distribution,
                                           final Parameter parameter) {

        if (parameter == distribution.getFieldParameter()) {
            return new GradientProvider() {
                @Override
                public int getDimension() {
                    return parameter.getDimension();
                }

                @Override
                public double[] getGradientLogDensity(Object x) {
                    return distribution.getGradientLogDensity(x);
                }
            };
        }

        final Parameter meanParameter = distribution.getMeanParameter();
        if (parameter == meanParameter) {
            final double[] scalarMeanGradient = meanParameter.getDimension() == 1 ? new double[1] : null;
            return new GradientProvider() {
                @Override
                public int getDimension() {
                    return meanParameter.getDimension();
                }

                @Override
                public double[] getGradientLogDensity(Object x) {
                    double[] precisionDiff = distribution.getPrecisionDiff((double[]) x);
                    if (meanParameter.getDimension() == distribution.getDimension()) {
                        return precisionDiff;
                    } else if (meanParameter.getDimension() == 1) {
                        double sum = 0.0;
                        for (int i = 0; i < distribution.getDimension(); ++i) {
                            sum += precisionDiff[i];
                        }
                        scalarMeanGradient[0] = sum;
                        return scalarMeanGradient;
                    } else {
                        throw new IllegalArgumentException("Unknown mean parameter structure");
                    }
                }
            };
        }

        if (parameter == distribution.getNuggetParameter()) {
            throw new RuntimeException("Not implemented");
        }

        for (BasisDimension basis : distribution.getBases()) {
            DesignMatrix designMatrix1 = basis.getDesignMatrix1();
            DesignMatrix designMatrix2 = basis.getDesignMatrix2();
            if (parameter == designMatrix1 || parameter == designMatrix2) {
                return new GradientProvider() {
                    @Override
                    public int getDimension() {
                        return parameter.getDimension();
                    }

                    @Override
                    public double[] getGradientLogDensity(Object x) {
                        throw new RuntimeException("Not yet implemented (DesignMatrix Gradient)");
                    }
                };
            }

            for (Parameter kernelParameter : basis.getKernel().getParameters()) {
                if (parameter == kernelParameter) {
                    throw new RuntimeException("Use GaussianProcessKernelGradient");
                }
            }
        }

        throw new IllegalArgumentException("Unknown parameter");
    }

    static GaussianProcessField getGaussianProcessField(RandomField randomField) {
        if (randomField.getDistribution() instanceof GaussianProcessField) {
            return (GaussianProcessField) randomField.getDistribution();
        }
        throw new IllegalArgumentException("Random field distribution is not a GaussianProcessField");
    }
}

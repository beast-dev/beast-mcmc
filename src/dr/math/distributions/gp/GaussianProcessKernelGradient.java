package dr.math.distributions.gp;

import dr.inference.distribution.RandomField;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.DesignMatrix;
import dr.inference.model.GradientProvider;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.xml.Reportable;
import java.util.List;

public class GaussianProcessKernelGradient implements Reportable, GradientWrtParameterProvider {

    private final Parameter field;
    private final int dim;
    private final AdditiveGaussianProcessDistribution distribution;
    private final GaussianProcessKernel.Base kernel;
    private final DesignMatrix designMatrix1;
    private final DesignMatrix designMatrix2;

//    private final boolean doScale;
//    private final boolean doLength;

    private final int parametersCount;
    private final GradientProvider provider;
    private double[] matrix;

    private final boolean DEBUG = true;

    public GaussianProcessKernelGradient(RandomField randomField,
                                         AdditiveGaussianProcessDistribution.BasisDimension basis,
                                         List<Parameter> parametersList) {
        this.field = randomField.getField();
        this.dim = field.getDimension();

        this.distribution = (AdditiveGaussianProcessDistribution) randomField.getDistribution();
        this.kernel = (GaussianProcessKernel.Base) basis.getKernel();

        // TODO allow for designMatrices to change?
        this.designMatrix1 = basis.getDesignMatrix1();
        this.designMatrix2 = basis.getDesignMatrix2();
        this.matrix = new double[dim * dim];

        List<Parameter> kernelParameters = kernel.getParameters();
        boolean doScale = false;
        boolean doLength = false;
        int k = 0;
        for (Parameter p : parametersList) {
            if (p == kernelParameters.get(0)) {
                doScale = true;
                k++;
            } else if (p == kernelParameters.get(1)) {
                doLength = true;
                k++;
            }
        }
        this.parametersCount = k;

        if (parametersCount == 1) {
            if (doScale) {
                this.provider = new GradientProvider() {
                    @Override
                    public int getDimension() {
                        return 1;
                    }

                    @Override
                    public double[] getGradientLogDensity(Object x) {
                        return new double[]{getGradientScale()};
                    }
                };
            } else if (doLength) {
                this.provider = new GradientProvider() {
                    @Override
                    public int getDimension() {
                        return 1;
                    }

                    @Override
                    public double[] getGradientLogDensity(Object x) {
                        return new double[]{getGradientLength()};
                    }
                };
            } else {
                throw new IllegalArgumentException("Gradient implemented only for scale or length");
            }
        } else if (parametersCount == 2) {
           this.provider = new GradientProvider() {
                @Override
                public int getDimension() {
                    return 2;
                }

                @Override
                public double[] getGradientLogDensity(Object x) {
                    double[] gradient = new double[2];
                    gradient[0] = getGradientScale();
                    gradient[1] = getGradientLength();
                    return gradient;
                }
            };
        } else {
            throw new IllegalArgumentException("Gradient not implemented for more than 2 parameters");
        }
    }

    public double getGradientScale() {
        if (DEBUG) System.out.println("Gradient Scale");
        return getGeneralGradient(kernel.getScale(), kernel.getScaleGradientFunction());
    }

    public double getGradientLength() {
        if (DEBUG) System.out.println("Gradient Length");
        return getGeneralGradient(kernel.getLength(), kernel.getLengthGradientFunction());
    }

    @Override
    public Likelihood getLikelihood() {
        return null;
    }

    @Override
    public Parameter getParameter() {
        return null;
    }

    @Override
    public int getDimension() {
        return parametersCount;
    }

    @Override
    public double[] getGradientLogDensity() {
        return provider.getGradientLogDensity(0.0);
    }


    private double getGeneralGradient(double hyperValue,
                                      GaussianProcessKernel.HyperparameterGradientFunction gradientFunction) {
        final double[] P = distribution.getPrecision();
        final double[] precisionDiff = distribution.getPrecisionDiff(field.getParameterValues());

        // Compute \frac{\partial K}{\partial \theta}
        for (int i = 0; i < dim; ++i) {
            double a = designMatrix1.getParameterValue(i, 0);
            for (int j = 0; j < dim; ++j) {
                double b = designMatrix2.getParameterValue(j, 0);
                matrix[i * dim + j] = gradientFunction.apply(a, b, hyperValue);
            }
        }
        return computeGeneralGradient(precisionDiff, P, matrix, dim);
    }

    private double computeGeneralGradient(double[] alpha, double[] P, double[] matrix, int dim) { // notice that the sign of alpha does not matter

        double quadForm = 0.0;
        double traceAB = 0.0;

        // Compute alpha^T B alpha
        // Compute trace(AB) = sum_{i,j} A_{ij} * B_{ij}
        int idx = 0;
        for (int i = 0; i < dim; ++i) {
            double alpha_i = alpha[i];
            for (int j = 0; j < dim; ++j, ++idx) {
                double alpha_j = alpha[j];
                double matrixij = matrix[idx];
                quadForm += alpha_i * matrixij * alpha_j;
                traceAB += P[idx] * matrixij;
            }
        }

        return 0.5 * (quadForm - traceAB);
    }

    @Override
    public String getReport() {
        StringBuilder sb = new StringBuilder();
        sb.append("gradient:");
        for (double value : getGradientLogDensity()) {
            sb.append(" ").append(value);
        }
        sb.append("\n");
        return sb.toString();
    }

}

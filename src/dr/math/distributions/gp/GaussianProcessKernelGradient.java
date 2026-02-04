package dr.math.distributions.gp;

import dr.inference.distribution.RandomField;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.*;
import dr.xml.Reportable;
import java.util.List;

public class GaussianProcessKernelGradient implements GradientWrtParameterProvider, Reportable {

    private final RandomField randomField;
    private final Parameter field;
    private final int dim;
    private final AdditiveGaussianProcessDistribution distribution;
    private final GaussianProcessKernel.Base kernel;
    private final DesignMatrix designMatrix1;
    private final DesignMatrix designMatrix2;

    private final GradientProvider provider;
    private final int parametersCount;
    private final Parameter hyperParameter;

    private final double[] matrix;

    private final boolean DEBUG = false;
    private  final BasisDimension basis;

    public GaussianProcessKernelGradient(RandomField randomField,
                                         BasisDimension basis,
                                         List<Parameter> parametersList) {
        this.randomField = randomField;
        this.field = randomField.getField();
        this.dim = field.getDimension();
        this.basis = basis;

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
            hyperParameter = parametersList.get(0);
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
            } else if (doLength) { // TODO This condition is always true
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
            hyperParameter = new CompoundParameter( null);
            for (int i = 0; i < parametersCount; i++) {
                ((CompoundParameter) hyperParameter).addParameter(kernelParameters.get(i));
            }
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
        return randomField;
    }

    @Override
    public Parameter getParameter() {
        return hyperParameter;
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
            double xi = designMatrix1.getParameterValue(i, 0);
            for (int j = 0; j < dim; ++j) {
                double xj = designMatrix2.getParameterValue(j, 0);
                double weight = 1.0;
                if (basis.getWeightFunction() != null) {
                    final double weightXi = basis.getWeightFunction().getWeight(xi);
                    final double weightXj = basis.getWeightFunction().getWeight(xj);
                    weight *= weightXi * weightXj;
                }
                matrix[i * dim + j] = gradientFunction.apply(xi, xj, hyperValue) * weight;
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
                double matrix_ij = matrix[idx];
                quadForm += alpha_i * matrix_ij * alpha_j;
                traceAB += P[idx] * matrix_ij;
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

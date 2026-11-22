package dr.math.distributions;

import dr.inference.model.Parameter;

public abstract class SplineGenerator {

    private final int degree;
    final int order;
    private final boolean intercept;
    private final boolean zeroOutOfBounds;

    public SplineGenerator(int degree, boolean intercept, boolean zeroOutOfBounds) {
        this.degree = degree;
        this.order = degree + 1;
        this.intercept = intercept;
        this.zeroOutOfBounds = zeroOutOfBounds;
    }

    public void fillBasis(double[][] basisMatrix,
                          double[] knots,
                          Parameter x,
                          double lower, double upper) {

        int offset = intercept ? 0 : 1;

        int rowDimension = basisMatrix.length;

        for (int r = 0; r < rowDimension; r++) {
            final double v = x.getParameterValue(r);
            if (Double.isFinite(v)) {
                if (v < lower || v > upper) {
                    if (zeroOutOfBounds) {
                        for (int i = offset; i < knots.length - order; ++i) {
                            basisMatrix[r][i - offset] = 0.0;
                        }
                    } else {
                        throw new RuntimeException("Out of spline bounds");
                    }
                } else {
                    for (int i = offset; i < knots.length - order; i++) {
                        basisMatrix[r][i - offset] = dispatch(i, degree, v, knots);
                    }
                }
            } else {
                for (int i = offset; i < knots.length - order; ++i) {
                    basisMatrix[r][i - offset] = Double.NaN;
                }
            }
        }
    }

    double getSplineBasis(int i, int k, double x, double[] knots) {

        if (k == 0) {
            if (knots[i] <= x && x < knots[i + 1]) {
                return 1.0;
            } else if (x == knots[knots.length - 1] && i == knots.length - order - 1) {
                return 1.0;
            } else {
                return 0.0;
            }
        }

        double denominator1 = knots[i + k] - knots[i];
        double denominator2 = knots[i + k + 1] - knots[i + 1];

        double term1 = 0.0;
        double term2 = 0.0;

        if (denominator1 != 0) {
            term1 = ((x - knots[i]) / denominator1) * getSplineBasis(i, k - 1, x, knots);
        }
        if (denominator2 != 0) {
            term2 = ((knots[i + k + 1] - x) / denominator2) * getSplineBasis(i + 1, k - 1, x, knots);
        }

        return term1 + term2;
    }

    abstract double dispatch(int i, int k, double x, double[] knots);

    static class Default extends SplineGenerator {

        public Default(int degree, boolean intercept, boolean zeroOutOfBounds) {
            super(degree, intercept, zeroOutOfBounds);
        }

        public double dispatch(int i, int k, double x, double[] knots) {
            return getSplineBasis(i, k, x, knots);
        }
    }

    public static class Derivative extends SplineGenerator {

        public Derivative(int degree, boolean intercept, boolean zeroOutOfBounds) {
            super(degree, intercept, zeroOutOfBounds);
        }

        public double dispatch(int i, int k, double x, double[] knots) {
            return getSplineBasisDerivative(i, k, x, knots);
        }

        double getSplineBasisDerivative(int i, int k, double x, double[] knots) {

            if (k == 0) {
                return 0.0;
            }

            double denominator1 = knots[i + k] - knots[i];
            double denominator2 = knots[i + k + 1] - knots[i + 1];

            double term1 = 0.0;
            double term2 = 0.0;

            if (denominator1 != 0) {
                term1 = getSplineBasis(i, k - 1, x, knots) / denominator1;
            }
            if (denominator2 != 0) {
                term2 = getSplineBasis(i + 1, k - 1, x, knots) / denominator2;
            }

            return k * (term1 - term2);
        }
    }
}

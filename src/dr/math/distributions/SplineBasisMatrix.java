package dr.math.distributions;

import dr.inference.model.*;
import dr.math.matrixAlgebra.WrappedMatrix;

import java.util.Arrays;

/**
 * A model that represents the B-spline basis matrix as a matrix parameter.
 * Refactored from BSplines.getBasisMatrixParameter().
 *
 * @author Pratyusa Datta
 * @author Marc Suchard
 */
public class SplineBasisMatrix extends MatrixParameter {

    private final Parameter x;
    private final Parameter k;
    private final int degree;
    private final int order;
    private final boolean intercept;

    private final Double lowerBoundary;
    private final Double upperBoundary;

    private double lower;
    private double upper;

    private final double[][] basisMatrix;
    private final double[] knots;

    private boolean basisMatrixKnown;
    private boolean knotsKnown;

    public SplineBasisMatrix(String name,
                             Parameter x,
                             Parameter k,
                             int degree,
                             boolean intercept,
                             Double lowerBoundary,
                             Double upperBoundary) {
        super(name, x.getDimension(), k.getDimension() + degree + (intercept ? 1 : 0));

        this.x = x;
        this.k = k;
        this.degree = degree;
        this.order = degree + 1;
        this.intercept = intercept;

        this.lowerBoundary = lowerBoundary;
        this.upperBoundary = upperBoundary;

        this.knots = new double[2 * order + k.getDimension()];
        this.basisMatrix = new double[rowDimension][columnDimension];

        x.addParameterListener(this);
        k.addParameterListener(this);

        basisMatrixKnown = false;
        knotsKnown = false;
    }

    private void getExpandedKnots() {

        if (!knotsKnown) {
            if (lowerBoundary != null) {
                lower = lowerBoundary;
            } else {
                lower = x.getParameterValue(0);
                for (int i = 1; i < x.getDimension(); ++i) {
                    double v = x.getParameterValue(i);
                    if (v < lower) {
                        lower = v;
                    }
                }
            }

            if (upperBoundary != null) {
                upper = upperBoundary;
            } else {
                upper = x.getParameterValue(0);
                for (int i = 0; i < x.getDimension(); ++i) {
                    double v = x.getParameterValue(i);
                    if (v > upper) {
                        upper = v;
                    }
                }
            }

            for (int i = 0; i < order; i++) {
                knots[i] = lower;
                knots[order + i] = upper;
            }

            for (int i = 0; i < k.getDimension(); ++i) {
                knots[2 * order + i] = k.getParameterValue(i);
            }

            Arrays.sort(knots);
            knotsKnown = true;
        }
    }

    private double getSplineBasis(int i, int d, double x, double[] knots) {
        if (d == 0) {
            if (knots[i] <= x && x < knots[i + 1]) {
                return 1.0;
            } else if (x == knots[knots.length - 1] && i == knots.length - order - 1) {
                return 1.0;
            } else {
                return 0.0;
            }
        }

        double denominator1 = knots[i + d] - knots[i];
        double denominator2 = knots[i + d + 1] - knots[i + 1];

        double term1 = 0.0;
        double term2 = 0.0;

        if (denominator1 != 0) {
            term1 = ((x - knots[i]) / denominator1) * getSplineBasis(i, d - 1, x, knots);
        }
        if (denominator2 != 0) {
            term2 = ((knots[i + d + 1] - x) / denominator2) * getSplineBasis(i + 1, d - 1, x, knots);
        }

        return term1 + term2;
    }

    private void computeBasisMatrix() {

        getExpandedKnots();

        if (!basisMatrixKnown) {

            int offset = intercept ? 0 : 1;

            for (int r = 0; r < rowDimension; r++) {
                final double v = x.getParameterValue(r);
                if (v < lower || v > upper) {
                    throw new RuntimeException("Out of spline bounds");
                }
                for (int i = offset; i < knots.length - order; i++) {
                    basisMatrix[r][i - offset] = getSplineBasis(i, degree, v, knots);
                }
            }
            basisMatrixKnown = true;
        }
    }

    @Override
    public double getParameterValue(int index) {
        int row = index / columnDimension;
        int col = index % columnDimension;

        return getParameterValue(row, col);
    }

    @Override
    public String getReport() {
        return new WrappedMatrix.MatrixParameter(this).toString();
    }

    @Override
    public double getParameterValue(int row, int col) {
        computeBasisMatrix();
        return basisMatrix[row][col];
    }

    @Override
    public void variableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        super.variableChangedEvent(variable, index, type);

        if (variable == x) {
            basisMatrixKnown = false;
        } else {
            throw new IllegalArgumentException("Unknown variable");
        }
    }

    @Override
    protected void storeValues() {
        super.storeValues();
    }

    @Override
    protected void restoreValues() {
        super.restoreValues();
    }
}
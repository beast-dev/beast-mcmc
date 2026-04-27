package dr.math.distributions;

import dr.inference.model.*;
import org.ejml.data.DenseMatrix64F;

import java.util.Arrays;

/**
 * A model that represents the B-spline basis matrix as a matrix parameter.
 * Refactored from BSplines.getBasisMatrixParameter().
 *
 * @author Pratyusa Datta
 * @author Marc Suchard
 */
public class SplineBasisMatrix extends DesignMatrix {

    private final Parameter x;
    private final Parameter k;
    private final int degree;
    private final int order;
    private final boolean intercept;
    private final boolean zeroOutOfBound;

    private final SplineGenerator generator;

    private final Double lowerBoundary;
    private final Double upperBoundary;

    private double lower;
    private double upper;

    private final double[][] basisMatrix;
    private final double[] knots;

    private boolean basisMatrixKnown;
    private boolean knotsKnown;
    private boolean variableChanged;

    public SplineBasisMatrix(String name,
                             Parameter x,
                             Parameter k,
                             int degree,
                             boolean intercept,
                             Double lowerBoundary,
                             Double upperBoundary,
                             boolean zeroOutOfBounds) {
        super(name, false);
        this.rowDimension = x.getDimension();
        this.columnDimension = (k != null ? k.getDimension() : 0) + degree + (intercept ? 1 : 0);

        this.x = x;
        this.k = k;
        this.degree = degree;
        this.order = degree + 1;
        this.intercept = intercept;
        this.zeroOutOfBound = zeroOutOfBounds;

        this.lowerBoundary = lowerBoundary;
        this.upperBoundary = upperBoundary;

        this.knots = new double[2 * order + (k != null ? k.getDimension() : 0)];
        this.basisMatrix = new double[rowDimension][columnDimension];

        this.generator = new SplineGenerator.Default(degree, intercept, zeroOutOfBounds);

        addParameter(x);

        if (k != null) {
            addParameter(k);
        }

        basisMatrixKnown = false;
        knotsKnown = false;
        variableChanged = false;
    }

    @Override
    public int getColumnDimension() {
        return columnDimension;
    }

    @Override
    public int getRowDimension() {
        return rowDimension;
    }

    private void computeExpandedKnots() {

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

            int kDim = (k != null ? k.getDimension() : 0);
            for (int i = 0; i < kDim; ++i) {
                knots[2 * order + i] = k.getParameterValue(i);
            }

            Arrays.sort(knots);
            knotsKnown = true;
        }
    }

    public Parameter getDesignParameter() {
        return x;
    }

    public Parameter getKnots() {
        return k;
    }

    public double[] getExpandedKnots() {
        computeExpandedKnots();
        return knots;
    }

    public double getLowerBoundary() {
        computeExpandedKnots();
        return lower;
    }

    public double getUpperBoundary() {
        computeExpandedKnots();
        return upper;
    }

    public int getDegree() { return degree; }

    public boolean getIncludeIntercept() { return intercept; }

    public boolean getZeroOutOfBound() { return zeroOutOfBound; }

    private void computeBasisMatrix() {

        computeExpandedKnots();

        if (!basisMatrixKnown) {
            generator.fillBasis(basisMatrix, knots, x, lower, upper);
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
        DenseMatrix64F mat = new DenseMatrix64F(getParameterAsMatrix());
        return mat.toString();
    }

    @Override
    public double getRawParameterValue(int row, int col) {
        computeBasisMatrix();
        return basisMatrix[row][col];
    }

    @Override
    public void variableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        super.variableChangedEvent(variable, index, type);

        if (variable == x) {
            basisMatrixKnown = false;
            variableChanged = true;
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
        if (variableChanged) {
            basisMatrixKnown = false;
            variableChanged = false;
        }
    }
}

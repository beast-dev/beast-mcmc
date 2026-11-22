package dr.math.distributions;

import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.loggers.NumberColumn;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.inference.model.VariableListener;

/**
 * A model that represents the B-spline basis matrix as a matrix parameter.
 * Refactored from BSplines.getBasisMatrixParameter().
 *
 * @author Pratyusa Datta
 * @author Marc Suchard
 */
public class SplinePredictionLogger implements Loggable, VariableListener {

    private final SplineBasisMatrix originalBasis;
    private final Parameter coefficients;
    private final Parameter predictionPoints;
    private final Parameter intercept;
    private final String name;
    private final String format;
    private final int dim;

    private SplineBasisMatrix predictionBasis;
    private LogColumn[] columns;

    public SplinePredictionLogger(String name,
                                  SplineBasisMatrix originalSplineBasis,
                                  Parameter coefficients,
                                  Parameter predictionPoints,
                                  Parameter intercept,
                                  String format) {

        this.name = name;
        this.originalBasis = originalSplineBasis;
        this.coefficients = coefficients;
        this.predictionPoints = predictionPoints;
        this.intercept = intercept;
        this.format = format;
        this.dim = predictionPoints.getDimension();

        if (originalSplineBasis.getColumnDimension() != coefficients.getDimension()) {
            throw new IllegalArgumentException("Basis dimension != coefficient dimension");
        }

        predictionPoints.addParameterListener(this);
    }

    @Override
    public LogColumn[] getColumns() {
        if (columns == null) {
            columns = new LogColumn[dim];
            for (int i = 0; i < dim; ++i) {
                String label = (format == null) ?
                        name + i :
                        name + "(" + String.format(format, predictionPoints.getParameterValue(i)) + ")";

                columns[i] = new InnerProductColumn(label, i);
            }
        }
        return columns;
    }

    class InnerProductColumn extends NumberColumn {

        private final int index;

        public InnerProductColumn(String label, int index) {
            super(label);
            this.index = index;
        }

        @Override
        public double getDoubleValue() {
            if (predictionBasis == null) {
                predictionBasis = makePredictionBasis(predictionPoints);
            }

            double innerProduct = 0;
            for (int j = 0; j < coefficients.getDimension(); ++j) {
                innerProduct += predictionBasis.getParameterValue(index, j) * coefficients.getParameterValue(j);
            }

            if (intercept != null) {
                innerProduct += intercept.getParameterValue(0);
            }
            
            return innerProduct;
        }
    }

    @Override
    public void variableChangedEvent(Variable variable, int index, Variable.ChangeType type) {
        if (variable == predictionPoints) {
            throw new IllegalArgumentException("Prediction points may not yet change");
        } else {
            throw new IllegalArgumentException("Unknown variable");
        }
    }

    private SplineBasisMatrix makePredictionBasis(Parameter x) {
        return new SplineBasisMatrix(name, x,
                originalBasis.getKnots(), originalBasis.getDegree(), originalBasis.getIncludeIntercept(),
                null, null, false);
    }
}

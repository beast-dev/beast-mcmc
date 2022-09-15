package dr.evomodel.treedatalikelihood.preorder;

import dr.evomodel.treedatalikelihood.continuous.cdi.PrecisionType;
import dr.math.matrixAlgebra.ReadableMatrix;
import dr.math.matrixAlgebra.ReadableVector;
import dr.math.matrixAlgebra.WrappedMatrix;
import dr.math.matrixAlgebra.WrappedVector;
import org.ejml.data.DenseMatrix64F;

/**
 * @author Marc A. Suchard
 */
public class WrappedNormalSufficientStatistics {

    private final WrappedVector mean;
    private final WrappedMatrix precision;
    private final WrappedMatrix variance;
    private final double precisionScalar;

    public WrappedNormalSufficientStatistics(WrappedVector mean,
                                             WrappedMatrix precision,
                                             WrappedMatrix variance) {
        this.mean = mean;
        this.precision = precision;
        this.variance = variance;
        this.precisionScalar = 1.0;
    }

    public WrappedNormalSufficientStatistics(double[] buffer,
                                             int index,
                                             int dim,
                                             DenseMatrix64F Pd,
                                             PrecisionType precisionType) {

        int partialOffset = (precisionType.getPartialsDimension(dim)) * index;
        this.mean = new WrappedVector.Raw(buffer, partialOffset, dim);
        if (precisionType == PrecisionType.SCALAR) {
            this.precision = new WrappedMatrix.Raw(Pd.getData(), 0, dim, dim);
            this.precisionScalar = buffer[partialOffset + dim];
            this.variance = null;
        } else {
            this.precisionScalar = 1.0;
            this.precision = new WrappedMatrix.Raw(buffer, partialOffset + dim, dim, dim);
            this.variance = new WrappedMatrix.Raw(buffer, partialOffset + dim + dim * dim, dim, dim);
        }
    }

    public WrappedVector getMean() { return mean; }

    public WrappedMatrix getPrecision() { return precision; }

    public WrappedMatrix getVariance() { return variance; }

    public double getMean(int row) {
        return mean.get(row);
    }

    public double getPrecision(int row, int col) {
        return precision.get(row, col);
    }

    public double getPrecisionScalar() { return precisionScalar; }

    public String toString() {
        return mean + " " + precision;
    }
}

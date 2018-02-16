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
public class WrappedMeanPrecision {

    private final WrappedVector mean;
    private final WrappedMatrix precision;
    private final double precisionScalar;

    WrappedMeanPrecision(double[] buffer,
                         int index,
                         int dim,
                         DenseMatrix64F Pd,
                         PrecisionType precisionType) {

        int partialOffset = (dim + precisionType.getMatrixLength(dim)) * index;
        this.mean = new WrappedVector.Raw(buffer, partialOffset, dim);
        if (precisionType == PrecisionType.SCALAR) {
            this.precision = new WrappedMatrix.Raw(Pd.getData(), 0, dim, dim);
            this.precisionScalar = buffer[partialOffset + dim];
        } else {
            throw new RuntimeException("Not yet implemented");
        }
    }

    public ReadableVector getMean() { return mean; }

    public ReadableMatrix getPrecision() { return precision; }

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

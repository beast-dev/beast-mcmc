package dr.evomodel.treedatalikelihood.preorder;

import dr.evomodel.treedatalikelihood.continuous.cdi.PrecisionType;
import dr.math.matrixAlgebra.missingData.MissingOps;
import org.ejml.data.DenseMatrix64F;

/**
 * @author Marc A. Suchard
 */
public class NormalSufficientStatistics {

    private final DenseMatrix64F mean;
    private final DenseMatrix64F precision;

    NormalSufficientStatistics(double[] buffer,
                                      int partialOffset,
                                      int dim,
                                      DenseMatrix64F Pd,
                                      PrecisionType precisionType) {

        this.mean = MissingOps.wrap(buffer, partialOffset, dim, 1);
        this.precision = DenseMatrix64F.wrap(dim, dim,
                precisionType.getScaledPrecision(buffer, partialOffset, Pd.data, dim));

    }

    @SuppressWarnings("unused")
    public NormalSufficientStatistics(DenseMatrix64F mean,
                                      DenseMatrix64F precision) {
        this.mean = mean;
        this.precision = precision;
    }

    public double getMean(int row) {
        return mean.get(row);
    }

    public double getPrecision(int row, int col) {
        return precision.unsafe_get(row, col);
    }

    public String toString() {
        return mean + " " + precision;
    }

    public String toVectorizedString() {
        return toVectorizedString(mean.getData()) + " " + toVectorizedString(precision.getData());
    }

    private static String toVectorizedString(double[] vector) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < vector.length - 1; ++i) {
            sb.append(vector[i]).append(" ");
        }
        sb.append(vector[vector.length - 1]);
        return sb.toString();
    }

}

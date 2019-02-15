package dr.math.matrixAlgebra.missingData;

import org.ejml.data.DenseMatrix64F;

public class PermutationIndices {

    final DenseMatrix64F matrix;
    final int dim;

    int zeroCount;
    int nonZeroFiniteCount;
    int infiniteCount;

    int[] nonZeroFiniteIndices;
    int[] zeroIndices;

    public PermutationIndices(DenseMatrix64F matrix) {

        this.matrix = matrix;
        dim = matrix.getNumCols();
        assert (dim == matrix.getNumRows());

        for (int i = 0; i < dim; ++i) {
            double diagonal = matrix.get(i, i);
            if (Double.isInfinite(diagonal)) {
                ++infiniteCount;
            } else if (diagonal == 0.0) {
                ++zeroCount;
            } else {
                ++nonZeroFiniteCount;
            }
        }
    }

    public int getNumberOfZeroDiagonals() {
        return zeroCount;
    }

    public int getNumberOfNonZeroFiniteDiagonals() {
        return nonZeroFiniteCount;
    }

    public int getNumberOfInfiniteDiagonals() {
        return infiniteCount;
    }

    public int[] getNonZeroFiniteIndices() {
        if (nonZeroFiniteIndices == null) {
            makeIndices();
        }
        return nonZeroFiniteIndices;
    }

    public int[] getZeroIndices() {
        if (zeroIndices == null) {
            makeIndices();
        }
        return zeroIndices;
    }

    private void makeIndices() {
        nonZeroFiniteIndices = new int[nonZeroFiniteCount];
        zeroIndices = new int[zeroCount];

        int zeroIndex = 0;
        int nonZeroFiniteIndex = 0;

        for (int i = 0; i < dim; ++i) {
            double diagonal = matrix.get(i, i);
            if (Double.isInfinite(diagonal)) {
                // Do nothing, unless we ever need infiniteIndices
            } else if (diagonal == 0.0) {
                zeroIndices[zeroIndex] = i;
                ++zeroIndex;
            } else {
                nonZeroFiniteIndices[nonZeroFiniteIndex] = i;
                ++nonZeroFiniteIndex;
            }
        }
    }
}

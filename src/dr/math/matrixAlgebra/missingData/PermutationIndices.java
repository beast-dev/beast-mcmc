package dr.math.matrixAlgebra.missingData;

import org.ejml.data.DenseMatrix64F;

public class PermutationIndices {

    private final DenseMatrix64F matrix;
    private final int dim;

    private int zeroCount;
    private int nonZeroFiniteCount;
    private int infiniteCount;

    private int[] nonZeroFiniteIndices;
    private int[] zeroIndices;
    private int[] infiniteIndices;

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

    @SuppressWarnings("unused")
    public int getNumberOfZeroDiagonals() {
        return zeroCount;
    }

    public int getNumberOfNonZeroFiniteDiagonals() {
        return nonZeroFiniteCount;
    }

    @SuppressWarnings("unused")
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

    @SuppressWarnings("unused")
    public int[] getInfiniteIndices() {
        if (infiniteIndices == null) {
            makeIndices();
        }
        return infiniteIndices;
    }

    private void makeIndices() {
        nonZeroFiniteIndices = new int[nonZeroFiniteCount];
        zeroIndices = new int[zeroCount];
        infiniteIndices = new int[infiniteCount];

        int zeroIndex = 0;
        int nonZeroFiniteIndex = 0;
        int infiniteIndex = 0;

        for (int i = 0; i < dim; ++i) {
            double diagonal = matrix.get(i, i);
            if (Double.isInfinite(diagonal)) {
                infiniteIndices[infiniteIndex] = i;
                ++infiniteIndex;
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

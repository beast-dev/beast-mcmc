package dr.matrix;

public class SparseSquareUpperTriangular extends SparseCompressedMatrix {

    public SparseSquareUpperTriangular(double[][] U) {
        super(U);

        checkStructure();
    }
    public SparseSquareUpperTriangular(int[] majorStarts, int[] minorIndices, double[] values,
                                 int nMajor) {
        super(majorStarts, minorIndices, values, nMajor, nMajor);

        checkStructure();
    }

    private void checkStructure() throws IllegalArgumentException {
        for (int i = 0; i < nMajor; ++i) {
            if (minorIndices[majorStarts[i]] < i) {
                throw new IllegalArgumentException("Not upper triangular");
            } else if (minorIndices[majorStarts[i]] > i) {
                throw new IllegalArgumentException("Rank deficient");
            }
        }
    }

    @SuppressWarnings("unused")
    public double[] backSolveMatrixVector(double[] b, int bOffset) {
        double[] x = new double[nMajor];
        backSolveInPlaceMatrixVector(x, 0, b, bOffset);
        return x;
    }

    public void backSolveInPlaceMatrixVector(double[] x, int xOffset,
                                             double[] b, int bOffset) {

        for (int i = nMajor - 1; i >= 0; --i) {

            double sum = 0.0;
            for (int j = majorStarts[i] + 1; j < majorStarts[i + 1]; ++j) {
                sum += values[j] * x[minorIndices[j] + xOffset];
            }

            x[i  + xOffset] = (b[i + bOffset] - sum) / values[majorStarts[i]];
        }
    }
}

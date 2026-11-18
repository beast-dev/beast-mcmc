package dr.inference.model;

/**
 * Shared row-major indexing for strict triangular coordinate vectors.
 */
final class StrictTriangularMatrixIndexing {

    private StrictTriangularMatrixIndexing() { }

    static int getStrictTriangularDimension(final int dimension) {
        if (dimension < 1) {
            throw new IllegalArgumentException("dimension must be positive");
        }
        return dimension * (dimension - 1) / 2;
    }

    static int lowerRowOffset(final int row) {
        if (row < 0) {
            throw new IllegalArgumentException("row must be non-negative");
        }
        return row * (row - 1) / 2;
    }

    static int upperRowOffset(final int dimension, final int row) {
        if (dimension < 1) {
            throw new IllegalArgumentException("dimension must be positive");
        }
        if (row < 0 || row >= dimension) {
            throw new IllegalArgumentException("row out of bounds");
        }
        return row * (2 * dimension - row - 1) / 2;
    }

    static int lowerIndex(final int row, final int col) {
        if (row <= col) {
            throw new IllegalArgumentException("lowerIndex requires row > col");
        }
        return lowerRowOffset(row) + col;
    }

    static int upperIndex(final int dimension, final int row, final int col) {
        if (row >= col) {
            throw new IllegalArgumentException("upperIndex requires row < col");
        }
        return upperRowOffset(dimension, row) + (col - row - 1);
    }
}

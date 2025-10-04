package dr.matrix;

import dr.math.matrixAlgebra.ReadableMatrix;
import dr.math.matrixAlgebra.WrappedMatrix;

import java.util.ArrayList;
import java.util.List;

public class SparseCompressedMatrix {

    int[] majorStarts;
    int[] minorIndices;
    double[] values;

    int nMajor;
    int nMinor;

    @SuppressWarnings("unused")
    public SparseCompressedMatrix(double[][] matrix) {
       this(new WrappedMatrix.ArrayOfArray(matrix));
    }

    @SuppressWarnings("unused")
    public SparseCompressedMatrix(double[] matrix, int offset, int nMajor, int nMinor) {
        this(new WrappedMatrix.Raw(matrix, offset, nMajor, nMinor));
    }

    public int getEntryCount() {
        return values.length;
    }

    public SparseCompressedMatrix(int[] majorStarts, int[] minorIndices, double[] values,
                                  int nMajor, int nMinor) {
        this.majorStarts = majorStarts;
        this.minorIndices = minorIndices;
        this.values = values;

        this.nMajor = nMajor;
        this.nMinor = nMinor;
    }

    public SparseCompressedMatrix(ReadableMatrix matrix) {

        this.nMajor = matrix.getMajorDim();
        this.nMinor = matrix.getMinorDim();

        this.majorStarts = new int[nMajor + 1];
        List<Integer> columnIndices = new ArrayList<>();
        List<Double> values = new ArrayList<>();

        int offset = 0;
        for (int i = 0; i < nMajor; ++i) {
            majorStarts[i] = offset;

            for (int j = 0; j < nMinor; ++j) {
                double entry = matrix.get(i, j);
                if (entry != 0.0) {
                    columnIndices.add(j);
                    values.add(entry);
                    ++offset;
                }
            }
        }
        majorStarts[nMajor] = offset;

        this.minorIndices = columnIndices.stream()
                .mapToInt(Integer::intValue)
                .toArray();

        this.values = values.stream()
                .mapToDouble(Double::doubleValue)
                .toArray();
    }

    public SparseCompressedMatrix makeCopy() {
        return new SparseCompressedMatrix(
                majorStarts.clone(), minorIndices.clone(), values.clone(),
                nMajor, nMinor);
    }

    public void copyTo(SparseCompressedMatrix destination) {
        System.arraycopy(majorStarts, 0, destination.majorStarts, 0, majorStarts.length);
        System.arraycopy(minorIndices, 0, destination.minorIndices, 0, minorIndices.length);
        System.arraycopy(values, 0, destination.values, 0, values.length);
        destination.nMajor = nMajor;
        destination.nMinor = nMinor;
    }

    public void swap(SparseCompressedMatrix other) {
        int[] starts = other.majorStarts;
        other.majorStarts = this.majorStarts;
        this.majorStarts = starts;

        int[] indices = other.minorIndices;
        other.minorIndices = this.minorIndices;
        this.minorIndices = indices;

        double[] values = other.values;
        other.values = this.values;
        this.values = values;

        int major = other.nMajor;
        other.nMajor = this.nMajor;
        this.nMajor = major;

        int minor = other.nMinor;
        other.nMinor = this.nMinor;
        this.nMinor = minor;
    }

    public double[][] makeDense() {
        double[][] result = new double[nMajor][nMinor];

        for (int i = 0; i < nMajor; ++i) {
            int begin = majorStarts[i];
            int end = majorStarts[i + 1];
            for (int j = begin; j < end; ++j) {
                result[i][minorIndices[j]] = values[j];
            }
        }

        return result;
    }

    public void multiplyInPlaceMatrixVector(double[] vector, int vectorOffset,
                                            double[] result, int resultOffset) {

        assert vector.length - vectorOffset >= nMajor;
        assert result.length - resultOffset >= nMajor;

        for (int i = 0; i < nMajor; ++i) {
            final int begin = majorStarts[i];
            final int end = majorStarts[i + 1];

            double sum = 0.0;

            for (int k = begin; k < end; ++k) {
                final int j = minorIndices[k];
                final double value = values[k];

                sum += value * vector[vectorOffset + j];
            }

            result[resultOffset + i] = sum;
        }
    }

    public double[] multiplyMatrixVector(double[] vector, int vectorOffset) {

        double[] result = new double[nMajor];
        multiplyInPlaceMatrixVector(vector, vectorOffset, result, 0);

        return result;
    }

    public double multiplyVectorTransposeMatrixVector(double[] left, int leftOffset,
                                                      double[] right, int rightOffset) {

        assert left.length - leftOffset >= nMajor;
        assert right.length - rightOffset >= nMajor;

        double outer = 0.0;

        for (int i = 0; i < nMajor; ++i) {
            final int begin = majorStarts[i];
            final int end = majorStarts[i + 1];

            double inner = 0.0;

            for (int k = begin; k < end; ++k) {
                final int j = minorIndices[k];
                final double value = values[k];

                inner += value * right[rightOffset + j];
            }

            outer += left[leftOffset + i] * inner;
        }

        return outer;
    }

    public double multiplyVecDiffTransposeMatrixVecDiff(double[] left, int leftOffset,
                                                        double[] right, int rightOffset) {

        assert left.length - leftOffset >= nMajor;
        assert right.length - rightOffset >= nMajor;

        double outer = 0.0;

        for (int i = 0; i < nMajor; ++i) {
            final int begin = majorStarts[i];
            final int end = majorStarts[i + 1];

            double inner = 0.0;

            for (int k = begin; k < end; ++k) {
                final int j = minorIndices[k];
                final double value = values[k];

                inner += value * (left[leftOffset + j] - right[rightOffset + j]);
            }

            outer += (left[leftOffset + i] - right[rightOffset + i]) * inner;
        }

        return outer;
    }

    @SuppressWarnings("unused")
    public static void multiplyInPlaceMatrixVector(SparseCompressedMatrix matrix,
                                                   double[] vector, int vectorOffset,
                                                   double[] result, int resultOffset) {
        matrix.multiplyInPlaceMatrixVector(vector, vectorOffset, result, resultOffset);
    }

    @SuppressWarnings("unused")
    public static double[] multiplyMatrixVector(SparseCompressedMatrix matrix,
                                                double[] vector, int vectorOffset) {
        return matrix.multiplyMatrixVector(vector, vectorOffset);
    }

    @SuppressWarnings("unused")
    public static double multipleVectorTransposeMatrixVector(SparseCompressedMatrix matrix,
                                                             double[] left, int leftOffset,
                                                             double[] right, int rightOffset) {
        return matrix.multiplyVectorTransposeMatrixVector(left, leftOffset, right, rightOffset);
    }
}

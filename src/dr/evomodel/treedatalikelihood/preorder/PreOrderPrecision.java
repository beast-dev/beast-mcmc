package dr.evomodel.treedatalikelihood.preorder;

import dr.evomodel.treedatalikelihood.continuous.cdi.PrecisionType;
import dr.math.matrixAlgebra.CholeskyDecomposition;
import dr.math.matrixAlgebra.missingData.MissingOps;
import dr.matrix.SparseCompressedMatrix;
import dr.matrix.SparseSquareUpperTriangular;
import org.ejml.data.DenseMatrix64F;

import static dr.math.matrixAlgebra.missingData.MissingOps.safeInvert2;

public interface PreOrderPrecision {

    double getPrecision(int row, int col);

    double getVariance(int row, int col);

    DenseMatrix64F getRawPrecision();

    DenseMatrix64F getRawVariance();

    @SuppressWarnings("unused")
    double[][] getRawVarianceCholesky();

    DenseMatrix64F getRawPrecisionCopy();

    DenseMatrix64F getRawVarianceCopy();

    String toVectorizedString();

    default double getScalar() { return 1.0; }

    class Dense implements PreOrderPrecision {

        public final DenseMatrix64F precision;
        private final int dim;

        DenseMatrix64F variance = null;
        double[][] varianceCholesky = null;

        public Dense(double[] buffer,
                     int index,
                     int dim,
                     PreOrderPrecision unscaled,
                     PrecisionType precisionType) {
            this(DenseMatrix64F.wrap(dim, dim,
                    precisionType.getScaledPrecision(buffer,
                            (precisionType.getPartialsDimension(dim)) * index,
                            unscaled.getRawPrecision().data, dim)));
        }

        public Dense(double[] matrices,
                     int index,
                     int dim) {
            this(MissingOps.wrap(matrices, (dim * dim) * index, dim, dim));
        }

        public Dense(DenseMatrix64F precision) {
            this.precision = precision;
            this.dim = precision.numCols;
        }

        public Dense(DenseMatrix64F precision, DenseMatrix64F variance) {
            this.precision = precision;
            this.variance = variance;
            this.dim = precision.numCols;
        }

        @Override
        public String toVectorizedString() {
            String string = NormalSufficientStatistics.toVectorizedString(precision);
            if (variance != null) {
                string += " " + NormalSufficientStatistics.toVectorizedString(variance);
            }
            return string;
        }

        @Override
        public double getPrecision(int row, int col) {
            return precision.unsafe_get(row, col);
        }

        @Override
        public double getVariance(int row, int col) {
            computeVariance();
            return variance.unsafe_get(row, col);
        }

        @Override
        public DenseMatrix64F getRawPrecision() {
            return precision;
        }

        @Override
        public DenseMatrix64F getRawVariance() {
            computeVariance();
            return variance;
        }

        @Override
        public double[][] getRawVarianceCholesky() {
            computeCholesky();
            return varianceCholesky;
        }

        @Override
        public DenseMatrix64F getRawPrecisionCopy() {
            return precision.copy();
        }

        public DenseMatrix64F getRawVarianceCopy() {
            computeVariance();
            return variance.copy();
        }

        private void computeCholesky() {
            if (varianceCholesky == null) {
                varianceCholesky = CholeskyDecomposition.execute(getRawVariance().data, 0, dim);
            }
        }

        private void computeVariance() {
            if (variance == null) {
                variance = new DenseMatrix64F(precision.numRows, precision.numCols);
                safeInvert2(precision, variance, false);
            }
        }


    }

    class Sparse implements PreOrderPrecision {

        private final SparseCompressedMatrix matrix;
        private final SparseSquareUpperTriangular decomposition;
        private final double scalar;

        public Sparse(SparseCompressedMatrix matrix,
                      SparseSquareUpperTriangular decomposition,
                      double scalar) {
            this.matrix = matrix;
            this.decomposition = decomposition;
            this.scalar = scalar;
        }

        // TODO ABCD
        public Sparse(double[] buffer,
                      int index,
                      int dim,
                      PreOrderPrecision Pd,
                      PrecisionType precisionType) {
            this.matrix = ((PreOrderPrecision.Sparse) Pd).getPrecision();
            this.decomposition = ((PreOrderPrecision.Sparse) Pd).getCholeskyDecomposition();
            this.scalar = precisionType.getPrecisionScalar(buffer, index, dim);
        }

        @Override
        public double getPrecision(int row, int col) {
            throw new RuntimeException("Not implemented");
        }

        @Override
        public double getVariance(int row, int col) {
            throw new RuntimeException("Not implemented");
        }

        @Override
        public DenseMatrix64F getRawPrecision() {
            return null;
        }

        @Override
        public DenseMatrix64F getRawVariance() {
            return null;
        }

        @Override
        public double[][] getRawVarianceCholesky() {
            return null;
        }

        @Override
        public DenseMatrix64F getRawPrecisionCopy() {
            return null;
        }

        @Override
        public DenseMatrix64F getRawVarianceCopy() {
            return null;
        }

        @Override
        public String toVectorizedString() {
            return null;
        }

        public SparseCompressedMatrix getPrecision() {
            return matrix;
        }

        @Override
        public double getScalar() {
            return scalar;
        }

        public SparseSquareUpperTriangular getCholeskyDecomposition() {
            return decomposition;
        }
    }
}

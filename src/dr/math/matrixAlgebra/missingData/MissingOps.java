package dr.math.matrixAlgebra.missingData;

import dr.inference.model.MatrixParameterInterface;
import dr.math.matrixAlgebra.*;
import org.ejml.alg.dense.decomposition.lu.LUDecompositionAlt_D64;
import org.ejml.alg.dense.linsol.lu.LinearSolverLu_D64;
import org.ejml.alg.dense.misc.UnrolledDeterminantFromMinor;
import org.ejml.alg.dense.misc.UnrolledInverseFromMinor;
import org.ejml.data.DenseMatrix64F;
import org.ejml.factory.DecompositionFactory;
import org.ejml.factory.LinearSolverFactory;
import org.ejml.interfaces.decomposition.SingularValueDecomposition;
import org.ejml.interfaces.linsol.LinearSolver;
import org.ejml.ops.CommonOps;
import org.ejml.ops.SingularOps;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static dr.math.matrixAlgebra.missingData.InversionResult.Code.*;
import static dr.util.EuclideanToInfiniteNormUnitBallTransform.projection;

/**
 * @author Marc A. Suchard
 */
public class MissingOps {

    public static DenseMatrix64F wrap(final double[] source, final int offset,
                                      final int numRows, final int numCols) {
        double[] buffer = new double[numRows * numCols];
        return wrap(source, offset, numRows, numCols, buffer);
    }

    public static DenseMatrix64F wrap(final double[] source, final int offset,
                                      final int numRows, final int numCols,
                                      final double[] buffer) {
        System.arraycopy(source, offset, buffer, 0, numRows * numCols);
        return DenseMatrix64F.wrap(numRows, numCols, buffer);
    }

    public static DenseMatrix64F wrap(MatrixParameterInterface A) {
        return wrap(A.getParameterValues(), 0, A.getRowDimension(), A.getColumnDimension());
    }

    public static DenseMatrix64F wrapDiagonal(final double[] source, final int offset,
                                              final int dim) {
        double[] buffer = new double[dim * dim];
        return wrapDiagonal(source, offset, dim, buffer);
    }

    public static DenseMatrix64F wrapSpherical(final double[] source, final int offset,
                                               final int dim) {
        double[] buffer = new double[dim * dim];
        return wrapSpherical(source, offset, dim, buffer);
    }

    public static DenseMatrix64F wrapSpherical(final double[] source, final int offset,
                                               final int dim,
                                               final double[] buffer) {
        fillSpherical(source, offset, dim, buffer);
        DenseMatrix64F res = DenseMatrix64F.wrap(dim, dim, buffer);
        CommonOps.transpose(res); // Column major.
        return res;
    }

    private static void fillSpherical(final double[] source, final int offset,
                                      final int dim, final double[] buffer) {
        for (int i = 0; i < dim; i++) {
            System.arraycopy(source, offset + i * (dim - 1),
                    buffer, i * dim, dim - 1);
            buffer[(i + 1) * dim - 1] = projection(source, offset + i * (dim - 1), dim - 1);
        }
    }

    public static DenseMatrix64F wrapDiagonal(final double[] source, final int offset,
                                              final int dim,
                                              final double[] buffer) {
        for (int i = 0; i < dim; ++i) {
            buffer[i * dim + i] = source[i];
        }
        return DenseMatrix64F.wrap(dim, dim, buffer);
    }

    public static DenseMatrix64F copy(ReadableMatrix source) {
        final int len = source.getDim();
        double[] buffer = new double[len];
        for (int i = 0; i < len; ++i) {
            buffer[i] = source.get(i);
        }
        return DenseMatrix64F.wrap(source.getMinorDim(), source.getMajorDim(), buffer);
    }

    public static void copy(DenseMatrix64F source, WritableMatrix destination) {
        final int len = destination.getDim();
        for (int i = 0; i < len; ++i) {
            destination.set(i, source.get(i));
        }
    }

    public static void gatherRowsAndColumns(final DenseMatrix64F source, final DenseMatrix64F destination,
                                            final int[] rowIndices, final int[] colIndices) {

        final int rowLength = rowIndices.length;
        final int colLength = colIndices.length;
        final double[] out = destination.getData();

        int index = 0;
        for (int i = 0; i < rowLength; ++i) {
            final int rowIndex = rowIndices[i];
            for (int j = 0; j < colLength; ++j) {
                out[index] = source.unsafe_get(rowIndex, colIndices[j]);
                ++index;
            }
        }
    }

    public static void scatterRowsAndColumns(final DenseMatrix64F source, final DenseMatrix64F destination,
                                             final int[] rowIdices, final int[] colIndices, final boolean clear) {
        if (clear) {
            Arrays.fill(destination.getData(), 0.0);
        }

        final int rowLength = rowIdices.length;
        final int colLength = colIndices.length;
        final double[] in = source.getData();

        int index = 0;
        for (int i = 0; i < rowLength; ++i) {
            final int rowIndex = rowIdices[i];
            for (int j = 0; j < colLength; ++j) {
                destination.unsafe_set(rowIndex, colIndices[j], in[index]);
                ++index;
            }
        }
    }

    public static void unwrap(final DenseMatrix64F source, final double[] destination, final int offset) {
        System.arraycopy(source.getData(), 0, destination, offset, source.getNumElements());
    }

    public static boolean anyDiagonalInfinities(DenseMatrix64F source) {
        boolean anyInfinities = false;
        for (int i = 0; i < source.getNumCols() && !anyInfinities; ++i) {
            if (Double.isInfinite(source.unsafe_get(i, i))) {
                anyInfinities = true;
            }
        }
        return anyInfinities;
    }

    public static boolean allFiniteDiagonals(DenseMatrix64F source) {
        boolean allFinite = true;

        final int length = source.getNumCols();
        for (int i = 0; i < length; ++i) {
            allFinite &= !Double.isInfinite(source.unsafe_get(i, i));
        }
        return allFinite;
    }

    public static int countFiniteDiagonals(DenseMatrix64F source) {
        final int length = source.getNumCols();

        int count = 0;
        for (int i = 0; i < length; ++i) {
            final double d = source.unsafe_get(i, i);
            if (!Double.isInfinite(d)) {
                ++count;
            }
        }
        return count;
    }

    public static int countZeroDiagonals(DenseMatrix64F source) {
        final int length = source.getNumCols();

        int count = 0;
        for (int i = 0; i < length; ++i) {
            final double d = source.unsafe_get(i, i);
            if (d == 0.0) {
                ++count;
            }
        }
        return count;
    }

    public static boolean allZeroDiagonals(DenseMatrix64F source) {
        final int length = source.getNumCols();

        for (int i = 0; i < length; ++i) {
            if (source.unsafe_get(i, i) != 0.0) {
                return false;
            }
        }
        return true;
    }

    public static void getFiniteDiagonalIndices(final DenseMatrix64F source, final int[] indices) {
        final int length = source.getNumCols();

        int index = 0;
        for (int i = 0; i < length; ++i) {
            final double d = source.unsafe_get(i, i);
            if (!Double.isInfinite(d)) {
                indices[index] = i;
                ++index;
            }
        }
    }

    public static int countFiniteNonZeroDiagonals(ReadableMatrix source) {
        final int length = source.getMajorDim();

        int count = 0;
        for (int i = 0; i < length; ++i) {
            final double d = source.get(i, i);
            if (!Double.isInfinite(d) && d != 0.0) {
                ++count;
            }
        }
        return count;
    }

    public static int countFiniteNonZeroDiagonals(DenseMatrix64F source) {
        final int length = source.getNumCols();

        int count = 0;
        for (int i = 0; i < length; ++i) {
            final double d = source.unsafe_get(i, i);
            if (!Double.isInfinite(d) && d != 0.0) {
                ++count;
            }
        }
        return count;
    }

    public static void getFiniteNonZeroDiagonalIndices(final DenseMatrix64F source, final int[] indices) {
        final int length = source.getNumCols();

        int index = 0;
        for (int i = 0; i < length; ++i) {
            final double d = source.unsafe_get(i, i);
            if (!Double.isInfinite(d) && d != 0.0) {
                indices[index] = i;
                ++index;
            }
        }
    }

    public static void addToDiagonal(DenseMatrix64F source, double increment) {
        final int width = source.getNumRows();
        for (int i = 0; i < width; ++i) {
            source.unsafe_set(i, i, source.unsafe_get(i, i) + increment);
        }
    }

    public static double det(DenseMatrix64F mat) {
        int numCol = mat.getNumCols();
        int numRow = mat.getNumRows();
        if (numCol != numRow) {
            throw new IllegalArgumentException("Must be a square matrix.");
        } else if (numCol <= 6) {
            return numCol >= 2 ? UnrolledDeterminantFromMinor.det(mat) : mat.get(0);
        } else {
            LUDecompositionAlt_D64 alg = new LUDecompositionAlt_D64();
            if (alg.inputModified()) {
                mat = mat.copy();
            }

            return !alg.decompose(mat) ? 0.0D : alg.computeDeterminant().real;
        }
    }

    public static double invertAndGetDeterminant(DenseMatrix64F mat, DenseMatrix64F result) {

        final int numCol = mat.getNumCols();
        final int numRow = mat.getNumRows();
        if (numCol != numRow) {
            throw new IllegalArgumentException("Must be a square matrix.");
        }

        if (numCol <= 5) {

            if (numCol >= 2) {
                UnrolledInverseFromMinor.inv(mat, result);
            } else {
                result.set(0, 1.0D / mat.get(0));
            }

            return numCol >= 2 ?
                    UnrolledDeterminantFromMinor.det(mat) :
                    mat.get(0);
        } else {

            LUDecompositionAlt_D64 alg = new LUDecompositionAlt_D64();
            LinearSolverLu_D64 solver = new LinearSolverLu_D64(alg);
            if (solver.modifiesA()) {
                mat = mat.copy();
            }

            if (!solver.setA(mat)) {
                return Double.NaN;
            }

            solver.invert(result);

            return alg.computeDeterminant().real;

        }
    }

    public static InversionResult safeDeterminant(DenseMatrix64F source, boolean invert) {
        final int finiteCount = countFiniteNonZeroDiagonals(source);

        InversionResult result;

        if (finiteCount == 0) {
            result = new InversionResult(NOT_OBSERVED, 0, 0);
        } else {
//            LinearSolver<DenseMatrix64F> solver = LinearSolverFactory.pseudoInverse(true);
//            solver.setA(source);
//
//            SingularValueDecomposition<DenseMatrix64F> svd = solver.getDecomposition();
//            double[] values = svd.getSingularValues();
//
//            if (values == null) {
//                throw new RuntimeException("Unable to perform SVD");
//            }

            SingularValueDecomposition<DenseMatrix64F> svd = DecompositionFactory.svd(source.getNumRows(), source.getNumCols(), false, false, false);
            if (!svd.decompose(source)) throw new RuntimeException("SVD decomposition failed");
            double[] values = svd.getSingularValues();

            double eps = SingularOps.singularThreshold(svd);

            int dim = 0;
            double det = 1;
            for (int i = 0; i < values.length; ++i) {
                final double lambda = values[i];
                if (lambda > 100000 * eps) {
                    det *= lambda;
                    ++dim;
                }
            }

            if (!invert) {
                det = 1.0 / det;
            }

            result = new InversionResult(dim == source.getNumCols() ? FULLY_OBSERVED : PARTIALLY_OBSERVED, dim, det);
        }

        return result;
    }

    public static InversionResult safeSolve(DenseMatrix64F A,
                                            WrappedVector b,
                                            WrappedVector x,
                                            boolean getDeterminat) {
        final int dim = b.getDim();

        assert (A.getNumRows() == dim && A.getNumCols() == dim);

        final DenseMatrix64F B = wrap(b.getBuffer(), b.getOffset(), dim, 1);
        final DenseMatrix64F X = new DenseMatrix64F(dim, 1);

        InversionResult ir = safeSolve(A, B, X, getDeterminat);


        for (int row = 0; row < dim; ++row) {
            x.set(row, X.unsafe_get(row, 0));
        }

        return ir;
    }

    public static InversionResult safeSolve(DenseMatrix64F A, DenseMatrix64F B, DenseMatrix64F X, boolean getDeterminant) {

        final int finiteCount = countFiniteNonZeroDiagonals(A);

        InversionResult result;
        if (finiteCount == 0) {
            Arrays.fill(X.getData(), 0);
            result = new InversionResult(NOT_OBSERVED, 0, 0);
        } else {

            LinearSolver<DenseMatrix64F> solver = LinearSolverFactory.pseudoInverse(true);
            solver.setA(A);
            solver.solve(B, X);

            int dim = 0;
            double det = 1;

            if (getDeterminant) {
//                SingularValueDecomposition<DenseMatrix64F> svd = solver.getDecomposition();
//                double[] values = svd.getSingularValues();

                SingularValueDecomposition<DenseMatrix64F> svd = DecompositionFactory.svd(A.getNumRows(), A.getNumCols(), false, false, false);
                if (!svd.decompose(A)) throw new RuntimeException("SVD decomposition failed");
                double[] values = svd.getSingularValues();

                double eps = SingularOps.singularThreshold(svd);

                for (int i = 0; i < values.length; ++i) {
                    final double lambda = values[i];
                    if (lambda > 100000 * eps) {
                        det *= lambda;
                        ++dim;
                    }
                }
            }

            result = new InversionResult(dim == A.getNumCols() ? FULLY_OBSERVED : PARTIALLY_OBSERVED, dim, det);
        }

        return result;
    }

//    public static InversionResult safeInvert(ReadableMatrix source, WritableMatrix destination, boolean getDeterminant) {
//
//        final int dim = source.getMajorDim();
//        final int finiteCount = countFiniteNonZeroDiagonals(source);
//        double det = 0;
//
//        if (finiteCount == dim) {
//
//            DenseMatrix64F result = new DenseMatrix64F(dim, dim);
//            DenseMatrix64F copyOfSource = copy(source);
//            if (getDeterminant) {
//                det = invertAndGetDeterminant(copyOfSource, result);
//            } else {
//                CommonOps.invert(copyOfSource, result);
//            }
//
//            copy(result, destination);
//
//            return new InversionResult(FULLY_OBSERVED, dim, det);
//        }
//
//        return null;
//    }


//TODO: deprecate class in favor of lazy PermutationIndices

    private static class MissingPartition {
        final int[] fInds;
        final int[] zInds;
        final int[] infInds;

        private MissingPartition(DenseMatrix64F matrix) {
            int dim = matrix.numCols;
            assert (dim == matrix.numRows);
            final List<Integer> finiteIndices = new ArrayList<Integer>();
            final List<Integer> zeroIndices = new ArrayList<Integer>();
            final List<Integer> infiniteIndices = new ArrayList<Integer>();
            for (int i = 0; i < dim; i++) {
                double x = matrix.get(i, i);
                if (x == 0) {
                    zeroIndices.add(i);
                } else if (Double.isInfinite(x)) {
                    infiniteIndices.add(i);
                } else {
                    finiteIndices.add(i);
                }
            }

            this.fInds = new int[finiteIndices.size()];
            for (int i = 0; i < finiteIndices.size(); i++) {
                fInds[i] = finiteIndices.get(i);
            }

            this.zInds = new int[zeroIndices.size()];
            for (int i = 0; i < zeroIndices.size(); i++) {
                zInds[i] = zeroIndices.get(i);
            }

            this.infInds = new int[infiniteIndices.size()];
            for (int i = 0; i < infiniteIndices.size(); i++) {
                infInds[i] = infiniteIndices.get(i);
            }

        }
    }

    public static InversionResult safeInvert(DenseMatrix64F source, DenseMatrix64F destination, boolean getDeterminant) {

        final int dim = source.getNumCols();
        final int finiteCount = countFiniteNonZeroDiagonals(source);
        double det = 0;

        if (finiteCount == dim) {
            if (getDeterminant) {
                det = invertAndGetDeterminant(source, destination);
            } else {
                CommonOps.invert(source, destination);
            }
            return new InversionResult(FULLY_OBSERVED, dim, det);
        } else {
            if (finiteCount == 0) {
                Arrays.fill(destination.getData(), 0);
                return new InversionResult(NOT_OBSERVED, 0, 0);
            } else {
                final int[] finiteIndices = new int[finiteCount];
                getFiniteNonZeroDiagonalIndices(source, finiteIndices);

                final DenseMatrix64F subSource = new DenseMatrix64F(finiteCount, finiteCount);
                gatherRowsAndColumns(source, subSource, finiteIndices, finiteIndices);

                final DenseMatrix64F inverseSubSource = new DenseMatrix64F(finiteCount, finiteCount);
                if (getDeterminant) {
                    det = invertAndGetDeterminant(subSource, inverseSubSource);
                } else {
                    CommonOps.invert(subSource, inverseSubSource);
                }

                scatterRowsAndColumns(inverseSubSource, destination, finiteIndices, finiteIndices, true);

                return new InversionResult(PARTIALLY_OBSERVED, finiteCount, det);
            }
        }
    }

    //TODO: Just have one safeInvert function after checking to make sure it doesn't break anything
    // TODO: change all inversion to return logDeterminant
    public static InversionResult safeInvert2(DenseMatrix64F source, DenseMatrix64F destination, boolean getDeterminant) {

        final int dim = source.getNumCols();
        final PermutationIndices permutationIndices = new PermutationIndices(source);
        final int finiteCount = permutationIndices.getNumberOfNonZeroFiniteDiagonals();

        double det = 0;

        if (finiteCount == dim) {
            if (getDeterminant) {
                det = invertAndGetDeterminant(source, destination);
            } else {
                CommonOps.invert(source, destination);
            }
            return new InversionResult(FULLY_OBSERVED, dim, det);
        } else {
            if (finiteCount == 0) {
                Arrays.fill(destination.getData(), 0);
                return new InversionResult(NOT_OBSERVED, 0, 0);
            } else {

                final int[] finiteIndices = permutationIndices.getNonZeroFiniteIndices();
                final int[] zeroIndices = permutationIndices.getZeroIndices();

                final DenseMatrix64F subSource = new DenseMatrix64F(finiteCount, finiteCount);
                gatherRowsAndColumns(source, subSource, finiteIndices, finiteIndices);

                final DenseMatrix64F inverseSubSource = new DenseMatrix64F(finiteCount, finiteCount);
                if (getDeterminant) {
                    det = invertAndGetDeterminant(subSource, inverseSubSource);
                } else {
                    CommonOps.invert(subSource, inverseSubSource);
                }

                scatterRowsAndColumns(inverseSubSource, destination, finiteIndices, finiteIndices, true);

                for (int i = 0; i < zeroIndices.length; i++) {
                    int index = zeroIndices[i];
                    destination.set(index, index, Double.POSITIVE_INFINITY);
                }
                
                return new InversionResult(PARTIALLY_OBSERVED, finiteCount, det);
            }
        }
    }

    public static void matrixVectorMultiple(final DenseMatrix64F A,
                                            final WrappedVector x,
                                            final WrappedVector y,
                                            final int dim) {
        if (buffer.length < dim) {
            buffer = new double[dim];
        }

        for (int row = 0; row < dim; ++row) {
            double sum = 0.0;
            for (int col = 0; col < dim; ++col) {
                sum += A.unsafe_get(row, col) * x.get(col);
            }
            buffer[row] = sum;
        }

        for (int col = 0; col < dim; ++col) {
            y.set(col, buffer[col]);
        }
    }

    private static double[] buffer = new double[16];

    public static void safeWeightedAverage(final WrappedVector mi,
                                           final DenseMatrix64F Pi,
                                           final WrappedVector mj,
                                           final DenseMatrix64F Pj,
                                           final WrappedVector mk,
                                           final DenseMatrix64F Vk,
                                           final int dimTrait) {
//        countZeroDiagonals(Vk);
        final double[] tmp = new double[dimTrait];
        for (int g = 0; g < dimTrait; ++g) {
            double sum = 0.0;
            boolean iInf = Double.isInfinite(Pi.unsafe_get(g, g));
            boolean jInf = Double.isInfinite(Pj.unsafe_get(g, g));
            if (iInf && jInf) {
                throw new IllegalArgumentException("Both precision matrices are infinite in dimension " + g);
            } else if (iInf) {
                sum = mi.get(g);
            } else if (jInf) {
                sum = mj.get(g);
            } else {
                for (int h = 0; h < dimTrait; ++h) {
                    sum += Pi.unsafe_get(g, h) * mi.get(h);
                    sum += Pj.unsafe_get(g, h) * mj.get(h);
                }
            }

            tmp[g] = sum;
        }

        for (int g = 0; g < dimTrait; ++g) {
            double sum = 0.0;
            if (Vk.unsafe_get(g, g) == 0.0) {
                sum = tmp[g];
            } else {
                for (int h = 0; h < dimTrait; ++h) {
                    sum += Vk.unsafe_get(g, h) * tmp[h];
                }
            }
            mk.set(g, sum);
        }
    }

    public static void weightedAverage(final ReadableVector mi,
                                       final DenseMatrix64F Pi,
                                       final ReadableVector mj,
                                       final DenseMatrix64F Pj,
                                       final WritableVector mk,
                                       final DenseMatrix64F Vk,
                                       final int dimTrait) {
        final double[] tmp = new double[dimTrait];
        for (int g = 0; g < dimTrait; ++g) {
            double sum = 0.0;
            for (int h = 0; h < dimTrait; ++h) {
                sum += Pi.unsafe_get(g, h) * mi.get(h);
                sum += Pj.unsafe_get(g, h) * mj.get(h);
            }
            tmp[g] = sum;
        }
        for (int g = 0; g < dimTrait; ++g) {
            double sum = 0.0;
            for (int h = 0; h < dimTrait; ++h) {
                sum += Vk.unsafe_get(g, h) * tmp[h];
            }
            mk.set(g, sum);
        }
    }

    public static void weightedAverage(final ReadableVector mi,
                                       final ReadableMatrix Pi,
                                       final ReadableVector mj,
                                       final ReadableMatrix Pj,
                                       final WritableVector mk,
                                       final ReadableMatrix Vk,
                                       final int dimTrait) {
        final double[] tmp = new double[dimTrait];
        for (int g = 0; g < dimTrait; ++g) {
            double sum = 0.0;
            for (int h = 0; h < dimTrait; ++h) {
                sum += Pi.get(g, h) * mi.get(h);
                sum += Pj.get(g, h) * mj.get(h);
            }
            tmp[g] = sum;
        }
        for (int g = 0; g < dimTrait; ++g) {
            double sum = 0.0;
            for (int h = 0; h < dimTrait; ++h) {
                sum += Vk.get(g, h) * tmp[h];
            }
            mk.set(g, sum);
        }
    }

    public static void weightedAverage(final double[] ipartial,
                                       final int ibo,
                                       final DenseMatrix64F Pi,
                                       final double[] jpartial,
                                       final int jbo,
                                       final DenseMatrix64F Pj,
                                       final double[] kpartial,
                                       final int kbo,
                                       final DenseMatrix64F Vk,
                                       final int dimTrait) {
        final double[] tmp = new double[dimTrait];
        weightedAverage(ipartial, ibo, Pi, jpartial, jbo, Pj, kpartial, kbo, Vk, dimTrait, tmp);
    }

    public static void weightedSum(final double[] ipartial,
                                   final int ibo,
                                   final DenseMatrix64F Pi,
                                   final double[] jpartial,
                                   final int jbo,
                                   final DenseMatrix64F Pj,
                                   final int dimTrait,
                                   final double[] out) {
        for (int g = 0; g < dimTrait; ++g) {
            double sum = 0.0;
            for (int h = 0; h < dimTrait; ++h) {
                sum += Pi.unsafe_get(g, h) * ipartial[ibo + h];
                sum += Pj.unsafe_get(g, h) * jpartial[jbo + h];
            }
            out[g] = sum;
        }
    }

    public static void weightedSumActualized(final double[] ipartial,
                                             final int ibo,
                                             final DenseMatrix64F Pi,
                                             final double[] iactualization,
                                             final int ido,
                                             final double[] jpartial,
                                             final int jbo,
                                             final DenseMatrix64F Pj,
                                             final double[] jactualization,
                                             final int jdo,
                                             final int dimTrait,
                                             final double[] out) {
        for (int g = 0; g < dimTrait; ++g) {
            double sum = 0.0;
            for (int h = 0; h < dimTrait; ++h) {
                sum += iactualization[ido + g] * Pi.unsafe_get(g, h) * ipartial[ibo + h];
                sum += jactualization[jdo + g] * Pj.unsafe_get(g, h) * jpartial[jbo + h];
            }
            out[g] = sum;
        }
    }

    public static void weightedAverage(final double[] ipartial,
                                       final int ibo,
                                       final DenseMatrix64F Pi,
                                       final double[] jpartial,
                                       final int jbo,
                                       final DenseMatrix64F Pj,
                                       final double[] kpartial,
                                       final int kbo,
                                       final DenseMatrix64F Vk,
                                       final int dimTrait,
                                       final double[] tmp) {
        weightedSum(ipartial, ibo, Pi, jpartial, jbo, Pj, dimTrait, tmp);
        for (int g = 0; g < dimTrait; ++g) {
            double sum = 0.0;
            for (int h = 0; h < dimTrait; ++h) {
                // TODO Fix conditions
                if (!Double.isInfinite(Vk.unsafe_get(g, g)) && !Double.isInfinite(Vk.unsafe_get(h, h))) {
                    sum += Vk.unsafe_get(g, h) * tmp[h];
                }

            }
            kpartial[kbo + g] = sum;
        }
    }

    public static double weightedInnerProduct(final double[] partials,
                                              final int bo,
                                              final DenseMatrix64F P,
                                              final int dimTrait) {
        double SS = 0;

        // vector-matrix-vector
        for (int g = 0; g < dimTrait; ++g) {
            final double ig = partials[bo + g];
            for (int h = 0; h < dimTrait; ++h) {
                final double ih = partials[bo + h];
                SS += ig * P.unsafe_get(g, h) * ih;
            }
        }

        return SS;
    }

    public static double weightedInnerProductOfDifferences(final double[] source1,
                                                           final int source1Offset,
                                                           final double[] source2,
                                                           final int source2Offset,
                                                           final DenseMatrix64F P,
                                                           final int dimTrait) {
        double SS = 0;
        for (int g = 0; g < dimTrait; ++g) {
            final double gDifference = source1[source1Offset + g] - source2[source2Offset + g];

            for (int h = 0; h < dimTrait; ++h) {
                final double hDifference = source1[source1Offset + h] - source2[source2Offset + h];

                SS += gDifference * P.unsafe_get(g, h) * hDifference;
            }
        }

        return SS;
    }


    public static double weightedThreeInnerProduct(final double[] ipartials,
                                                   final int ibo,
                                                   final DenseMatrix64F Pip,
                                                   final double[] jpartials,
                                                   final int jbo,
                                                   final DenseMatrix64F Pjp,
                                                   final double[] kpartials,
                                                   final int kbo,
                                                   final DenseMatrix64F Pk,
                                                   final int dimTrait) {

        // TODO Is it better to split into 3 separate calls to weightedInnerProduct?

        double SSi = 0;
        double SSj = 0;
        double SSk = 0;


        // vector-matrix-vector TODO in parallel
        for (int g = 0; g < dimTrait; ++g) {
            final double ig = ipartials[ibo + g];
            final double jg = jpartials[jbo + g];
            final double kg = kpartials[kbo + g];

            for (int h = 0; h < dimTrait; ++h) {
                final double ih = ipartials[ibo + h];
                final double jh = jpartials[jbo + h];
                final double kh = kpartials[kbo + h];


                SSi += ig * Pip.unsafe_get(g, h) * ih;
                SSj += jg * Pjp.unsafe_get(g, h) * jh;
                SSk += kg * Pk.unsafe_get(g, h) * kh;

            }
        }

        return SSi + SSj - SSk;
    }


    public static void add(ReadableMatrix p1,
                           ReadableMatrix p2,
                           WritableMatrix p12) {

        assert (p1.getDim() == p2.getDim());
        assert (p1.getDim() == p12.getDim());

        final int dim = p12.getDim();

        for (int i = 0; i < dim; ++i) {
            p12.set(i, p1.get(i) + p2.get(i));
        }
    }
}

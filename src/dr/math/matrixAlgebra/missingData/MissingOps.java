package dr.math.matrixAlgebra.missingData;

import dr.inference.model.MatrixParameterInterface;
import dr.math.distributions.MultivariateNormalDistribution;
import dr.math.matrixAlgebra.*;
import org.ejml.alg.dense.decomposition.lu.LUDecompositionAlt_D64;
import org.ejml.alg.dense.linsol.lu.LinearSolverLu_D64;
import org.ejml.alg.dense.linsol.svd.SolvePseudoInverseSvd;
import org.ejml.alg.dense.misc.UnrolledDeterminantFromMinor;
import org.ejml.alg.dense.misc.UnrolledInverseFromMinor;
import org.ejml.data.DenseMatrix64F;
import org.ejml.factory.DecompositionFactory;
import org.ejml.factory.LinearSolverFactory;
import org.ejml.interfaces.decomposition.SingularValueDecomposition;
import org.ejml.interfaces.linsol.LinearSolver;
import org.ejml.ops.CommonOps;
import org.ejml.ops.SingularOps;

import java.util.Arrays;

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

    public static DenseMatrix64F wrapDiagonalInverse(final double[] source, final int offset,
                                                     final int dim) {
        double[] buffer = new double[dim * dim];
        return wrapDiagonalInverse(source, offset, dim, buffer);
    }

    public static DenseMatrix64F wrapDiagonalInverse(final double[] source, final int offset,
                                                     final int dim,
                                                     final double[] buffer) {
        for (int i = 0; i < dim; ++i) {
            buffer[i * dim + i] = 1 / source[i];
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

    public static void copyRowsAndColumns(final DenseMatrix64F source, final DenseMatrix64F destination,
                                          final int[] rowIndices, final int[] colIndices, final boolean clear) {
        if (clear) {
            Arrays.fill(destination.getData(), 0.0);
        }
        for (int row : rowIndices) {
            for (int col : colIndices) {
                destination.unsafe_set(row, col, source.unsafe_get(row, col));
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

    public static void unwrapIdentity(final double[] destination, final int offset, final int dim) {
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < i; j++) {
                destination[offset + i * dim + j] = 0.0;
            }
            destination[offset + i * dim + i] = 1.0;
            for (int j = i + 1; j < dim; j++) {
                destination[offset + i * dim + j] = 0.0;
            }
        }
    }

    public static void blockUnwrap(final DenseMatrix64F block, final double[] destination,
                                   final int destinationOffset,
                                   final int offsetRow, final int offsetCol,
                                   final int nCols) {
        for (int i = 0; i < block.getNumRows(); i++) { // Rows
            for (int j = 0; j < block.getNumCols(); j++) {
                destination[destinationOffset + (i + offsetRow) * nCols + j + offsetCol] = block.get(i, j);
            }
        }
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

    public static double invertAndGetDeterminant(DenseMatrix64F mat, DenseMatrix64F result, boolean log) {

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

            double det = numCol >= 2 ?
                    UnrolledDeterminantFromMinor.det(mat) :
                    mat.get(0);
            return log ? Math.log(det) : det;

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

            return log ? computeLogDeterminant(alg) : alg.computeDeterminant().real;

        }
    }

    private static double computeLogDeterminant(LUDecompositionAlt_D64 alg) {
        int n = alg.getLU().getNumCols();
        if (n != alg.getLU().getNumRows()) {
            throw new IllegalArgumentException("Must be a square matrix.");
        } else {
            double logDet = 0;
            double[] dataLU = alg.getLU().getData();
            for (int i = 0; i < n * n; i += n + 1) {
                logDet += Math.log(Math.abs(dataLU[i]));
            }

            return logDet;
        }
    }

    public static InversionResult safeDeterminant(DenseMatrix64F source, boolean invert) {
        final int finiteCount = countFiniteNonZeroDiagonals(source);

        InversionResult result;

        if (finiteCount == 0) {
            result = new InversionResult(NOT_OBSERVED, 0, Double.NEGATIVE_INFINITY, true);
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
            if (!svd.decompose(source)) {
                if (SingularOps.rank(svd) == 0)
                    return new InversionResult(NOT_OBSERVED, 0, Double.NEGATIVE_INFINITY, true);
                throw new RuntimeException("SVD decomposition failed");
            }
            double[] values = svd.getSingularValues();
            double tol = SingularOps.singularThreshold(svd);
//            double tol = 0.0;

            int dim = 0;
            double logDet = 0;
            for (int i = 0; i < values.length; i++) {
                final double lambda = values[i];
                if (lambda > tol) {
                    logDet += Math.log(lambda);
                    ++dim;
                }
            }

            if (invert) {
                logDet = -logDet;
            }

            result = new InversionResult(dim == source.getNumCols() ? FULLY_OBSERVED : PARTIALLY_OBSERVED, dim, logDet, true);
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

    public static InversionResult safeSolve(DenseMatrix64F A, DenseMatrix64F B, DenseMatrix64F X, boolean getLogDeterminant) {

        final int finiteCount = countFiniteNonZeroDiagonals(A);

        InversionResult result;
        if (finiteCount == 0) {
            Arrays.fill(X.getData(), 0);
            result = new InversionResult(NOT_OBSERVED, 0, Double.NEGATIVE_INFINITY, true);
        } else {

            LinearSolver<DenseMatrix64F> solver = LinearSolverFactory.pseudoInverse(true);
            ((SolvePseudoInverseSvd) solver).setThreshold(1e-8); // TODO No magic numbers, define as static final somewhere obvsious
            solver.setA(A);
            solver.solve(B, X);

            int dim = 0;
            double logDet = 0;

            //TODO: getLogDeterminant should never be used any more
            assert !getLogDeterminant;
//            if (getLogDeterminant) {
////                SingularValueDecomposition<DenseMatrix64F> svd = solver.getDecomposition();
////                double[] values = svd.getSingularValues();
//
//                SingularValueDecomposition<DenseMatrix64F> svd = DecompositionFactory.svd(A.getNumRows(), A.getNumCols(), false, false, false);
//                if (!svd.decompose(A)) {
//                    if (SingularOps.rank(svd) == 0)
//                        return new InversionResult(NOT_OBSERVED, 0, Double.NEGATIVE_INFINITY, true);
//                    throw new RuntimeException("SVD decomposition failed");
//                }
//                double[] values = svd.getSingularValues();
//
//                double tol = SingularOps.singularThreshold(svd);
//
//                for (int i = 0; i < values.length; ++i) {
//                    final double lambda = values[i];
//                    if (lambda > tol) {
//                        logDet += Math.log(lambda);
//                        ++dim;
//                    }
//                }
//            }

            result = new InversionResult(dim == A.getNumCols() ? FULLY_OBSERVED : PARTIALLY_OBSERVED, dim, logDet, true);
        }

        return result;
    }

    public static InversionResult safeInvertPrecision(DenseMatrix64F source, DenseMatrix64F destination,
                                                      boolean getLogDeterminant) {
        return safeInvert2(source, destination, getLogDeterminant, false);
    }

    public static InversionResult safeInvertVariance(DenseMatrix64F source, DenseMatrix64F destination,
                                                     boolean getLogDeterminant) {
        return safeInvert2(source, destination, getLogDeterminant, true);
    }

    public static InversionResult safeInvert2(DenseMatrix64F source, DenseMatrix64F destination,
                                              boolean getLogDeterminant) {
        return safeInvert2(source, destination, getLogDeterminant, true);
    }

    private static InversionResult safeInvert2(DenseMatrix64F source, DenseMatrix64F destination,
                                               boolean getLogDeterminant,
                                               boolean isInputVariance) {

        final int dim = source.getNumCols();
        final PermutationIndices permutationIndices = new PermutationIndices(source);
        final int finiteNonZeroCount = permutationIndices.getNumberOfNonZeroFiniteDiagonals();

        double logDet = 0;

        if (finiteNonZeroCount == dim) {
            if (getLogDeterminant) {
                logDet = invertAndGetDeterminant(source, destination, true);
            } else {
//                CommonOps.invert(source, destination);
                symmPosDefInvert(source, destination);
            }
            return new InversionResult(FULLY_OBSERVED, dim, logDet, true);
        } else {
            if (finiteNonZeroCount == 0) {
                Arrays.fill(destination.getData(), 0);
                //TODO: should NOT_OBSERVED vs FULLY_OBSERVED depend on whether this is a variance vs precision matrix?
                int infCount = permutationIndices.getNumberOfInfiniteDiagonals();

                if (infCount == dim) { //All infinity on diagonals of original matrix

                    return isInputVariance ?
                            new InversionResult(NOT_OBSERVED, 0, Double.NEGATIVE_INFINITY, true) :
                            new InversionResult(FULLY_OBSERVED, dim, Double.POSITIVE_INFINITY, true);

                } else {

                    int zeroCount = permutationIndices.getNumberOfZeroDiagonals();

                    if (zeroCount == dim) { //All zero on diagonals of original matrix
                        for (int i = 0; i < dim; i++) {
                            destination.set(i, i, Double.POSITIVE_INFINITY);
                        }
                        return isInputVariance ?
                                new InversionResult(FULLY_OBSERVED, dim, Double.POSITIVE_INFINITY, true) :
                                new InversionResult(NOT_OBSERVED, 0, Double.NEGATIVE_INFINITY, true);

                    } else { //Both zeros and infinities (but no non-zero finite entries) on diagonal
                        int[] zeroInds = permutationIndices.getZeroIndices();
                        for (int i : zeroInds) {
                            destination.set(i, i, Double.POSITIVE_INFINITY);
                        }

                        return isInputVariance ?
                                new InversionResult(PARTIALLY_OBSERVED, zeroCount, Double.POSITIVE_INFINITY, true) :
                                new InversionResult(PARTIALLY_OBSERVED, infCount, Double.POSITIVE_INFINITY, true);
                    }
                }

            } else {

                final int[] finiteIndices = permutationIndices.getNonZeroFiniteIndices();
                final int[] zeroIndices = permutationIndices.getZeroIndices();

                final DenseMatrix64F subSource = new DenseMatrix64F(finiteNonZeroCount, finiteNonZeroCount);
                gatherRowsAndColumns(source, subSource, finiteIndices, finiteIndices);

                final DenseMatrix64F inverseSubSource = new DenseMatrix64F(finiteNonZeroCount, finiteNonZeroCount);
                if (getLogDeterminant) {
                    logDet = invertAndGetDeterminant(subSource, inverseSubSource, true);
                } else {
//                    CommonOps.invert(subSource, inverseSubSource);
                    symmPosDefInvert(subSource, inverseSubSource);
                }

                scatterRowsAndColumns(inverseSubSource, destination, finiteIndices, finiteIndices, true);

                for (int i = 0; i < zeroIndices.length; i++) {
                    int index = zeroIndices[i];
                    destination.set(index, index, Double.POSITIVE_INFINITY);
                }

                int fullyObsCount = isInputVariance ? permutationIndices.getNumberOfZeroDiagonals() :
                        permutationIndices.getNumberOfInfiniteDiagonals();
                logDet = fullyObsCount == 0 ? logDet : Double.POSITIVE_INFINITY;


                return new InversionResult(PARTIALLY_OBSERVED, fullyObsCount + finiteNonZeroCount, logDet, true);
            }
        }
    }

    public static void symmPosDefInvert(DenseMatrix64F P, DenseMatrix64F P_inv) {
        LinearSolver<DenseMatrix64F> solver = LinearSolverFactory.symmPosDef(P.getNumCols());
        DenseMatrix64F Pbis = new DenseMatrix64F(P);
        if (!solver.setA(Pbis)) {
            CommonOps.invert(P, P_inv);
        } else {
            solver.invert(P_inv);
        }
    }

    public static void safeMult(DenseMatrix64F sourceA, DenseMatrix64F sourceB, DenseMatrix64F destination) {

        final int dim = sourceA.getNumCols();
        assert ((dim == sourceA.getNumRows()) && dim == sourceB.getNumCols() && dim == sourceB.getNumRows()) :
                "In safeMult, A and B must be square with the same dimension.";
        final PermutationIndices permutationIndicesA = new PermutationIndices(sourceA);
        final int infiniteCountA = permutationIndicesA.getNumberOfInfiniteDiagonals();
        final PermutationIndices permutationIndicesB = new PermutationIndices(sourceB);
        final int infiniteCountB = permutationIndicesB.getNumberOfInfiniteDiagonals();

        if (infiniteCountA == 0 && infiniteCountB == 0) {
            CommonOps.mult(sourceA, sourceB, destination);
        } else if (infiniteCountA == dim) {
            CommonOps.scale(1.0, sourceA, destination);
        } else if (infiniteCountB == dim) {
            CommonOps.scale(1.0, sourceB, destination);
        } else {
            throw new RuntimeException("Partial safeMult not yet implemented.");
        }
    }


//    public static void safeAdd(DenseMatrix64F source0, DenseMatrix64F source1, DenseMatrix64F destination) {
//        CommonOps.add(source0, source1, destination);
//
//        for (int i = 0; i < destination.numCols; ++i) {
//            if (Double.isInfinite(destination.unsafe_get(i, i))) {
//                for (int j = 0; j < destination.numRows; ++j) {
//                    if (i != j) {
//                        destination.unsafe_set(i, j, 0.0);
//                        destination.unsafe_set(j, i, 0.0);
//                    }
//                }
//            }
//        }
//    }

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
            if (!Double.isInfinite(Vk.unsafe_get(g, g))) {
                double sum = 0.0;
                for (int h = 0; h < dimTrait; ++h) {
                    if (!Double.isInfinite(Vk.unsafe_get(h, h))) {

                        sum += Vk.unsafe_get(g, h) * tmp[h];

                    }
                }
                kpartial[kbo + g] = sum;

            }
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

    public static double weightedThreeInnerProductNormalized(final double[] ipartials,
                                                             final int ibo,
                                                             final DenseMatrix64F Pip,
                                                             final double[] jpartials,
                                                             final int jbo,
                                                             final DenseMatrix64F Pjp,
                                                             final double[] kpartials,
                                                             final int kbo,
                                                             final double[] kpartialsBis,
                                                             final int kboBis,
                                                             final int dimTrait) {

        double SSi = 0;
        double SSj = 0;
        double SSk = 0;

        // vector-matrix-vector TODO in parallel
        for (int g = 0; g < dimTrait; ++g) {
            final double ig = ipartials[ibo + g];
            final double jg = jpartials[jbo + g];
            final double kg = kpartials[kbo + g];
            final double kgBis = kpartialsBis[kboBis + g];

            for (int h = 0; h < dimTrait; ++h) {
                final double ih = ipartials[ibo + h];
                final double jh = jpartials[jbo + h];

                SSi += ig * Pip.unsafe_get(g, h) * ih;
                SSj += jg * Pjp.unsafe_get(g, h) * jh;
            }
            SSk += kg * kgBis;
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

    public static void forceSymmetric(DenseMatrix64F P) {
        DenseMatrix64F Ptrans = new DenseMatrix64F(P);
        CommonOps.transpose(P, Ptrans);
        CommonOps.addEquals(P, Ptrans);
        CommonOps.scale(0.5, P);
    }

    public static void symmetricMult(DenseMatrix64F Q, DenseMatrix64F P, DenseMatrix64F QtPQ) {
        int dimTrait = Q.getNumCols();
        assert dimTrait == Q.getNumRows() && dimTrait == P.getNumCols() && dimTrait == P.getNumRows();
        for (int i = 0; i < dimTrait; i++) {
            for (int j = i; j < dimTrait; j++) {
                double val = 0;
                for (int k = 0; k < dimTrait; k++) {
                    for (int r = 0; r < dimTrait; r++) {
                        val += P.unsafe_get(k, r) * Q.unsafe_get(k, i) * Q.unsafe_get(r, j);
                    }
                }
                QtPQ.unsafe_set(i, j, val);
                QtPQ.unsafe_set(j, i, val);
            }
        }
    }

    public static void diagMult(double[] d, DenseMatrix64F M) {
        diagMult(d, M, M);
    }

    public static void diagMult(double[] d, DenseMatrix64F source, DenseMatrix64F dest) {
        assert d.length == source.getNumRows();
        assert source.getNumRows() == dest.getNumRows() && source.getNumCols() == dest.getNumCols();
        for (int i = 0; i < source.getNumRows(); i++) {
            for (int j = 0; j < source.getNumCols(); j++) {
                dest.unsafe_set(i, j, d[i] * source.unsafe_get(i, j));
            }
        }
    }

    public static void diagMult(DenseMatrix64F M, double[] d) {
        diagMult(M, d, M);
    }

    public static void diagMult(DenseMatrix64F source, double[] d, DenseMatrix64F dest) {
        assert d.length == source.getNumCols();
        assert source.getNumRows() == dest.getNumRows() && source.getNumCols() == dest.getNumCols();
        for (int i = 0; i < source.getNumRows(); i++) {
            for (int j = 0; j < source.getNumCols(); j++) {
                dest.unsafe_set(i, j, d[j] * source.unsafe_get(i, j));
            }
        }
    }

    public static void diagDiv(double[] d, DenseMatrix64F M) {
        assert d.length == M.getNumRows();
        for (int i = 0; i < M.getNumRows(); i++) {
            for (int j = 0; j < M.getNumCols(); j++) {
                M.unsafe_set(i, j, M.unsafe_get(i, j) / d[i]);
            }
        }
    }

    public static void diagDiv(DenseMatrix64F M, double[] d) {
        assert d.length == M.getNumCols();
        for (int i = 0; i < M.getNumRows(); i++) {
            for (int j = 0; j < M.getNumCols(); j++) {
                M.unsafe_set(i, j, M.unsafe_get(i, j) / d[j]);
            }
        }
    }

    public static void addTransEquals(DenseMatrix64F M) {
        assert M.getNumCols() == M.getNumRows();
        for (int i = 0; i < M.getNumCols(); i++) {
            M.unsafe_set(i, i, 2 * M.unsafe_get(i, i));
            for (int j = 0; j < i; j++) {
                M.unsafe_set(i, j, M.unsafe_get(i, j) + M.unsafe_get(j, i));
                M.unsafe_set(j, i, M.unsafe_get(i, j));
            }
        }
    }

    public static double[] nextPossiblyDegenerateNormal(ReadableVector mean, DenseMatrix64F variance) {
        int dim = mean.getDim();

        if (variance.numCols != dim || variance.numRows != dim) {
            throw new RuntimeException("Variance is a " + variance.numRows + "x" + variance.numCols +
                    " matrix but mean has dimension " + dim);
        }

        int zeroCount = countZeroDiagonals(variance);
        int nonZeroCount = countFiniteNonZeroDiagonals(variance);
        if (zeroCount + nonZeroCount != dim) {
            throw new RuntimeException("At least one diagonal element of the variance is infinity. " +
                    "Cannot sample from distribution with infinite variance");
        }


        double[] buffer = ReadableVector.Utils.toArray(mean);


        if (nonZeroCount == dim) {
            double[][] cholesky = CholeskyDecomposition.execute(variance.data, 0, dim);
            return MultivariateNormalDistribution.nextMultivariateNormalCholesky(buffer, cholesky);
        }

        int[] latentIndices = new int[nonZeroCount];

        int latI = 0;
        for (int i = 0; i < dim; i++) {
            if (variance.get(i, i) > 0) {
                latentIndices[latI] = i;
                latI++;
            }
        }

        WrappedMatrix.Indexed subVar = new WrappedMatrix.Indexed(variance.data, 0,
                latentIndices, latentIndices,
                dim, dim);


        WrappedVector.Indexed subMean = new WrappedVector.Indexed(buffer, 0, latentIndices);


        double[] latentDraw = MultivariateNormalDistribution.nextMultivariateNormalVariance(
                ReadableVector.Utils.toArray(subMean), ReadableMatrix.Utils.toMatrixArray(subVar));

        for (int i = 0; i < latentIndices.length; i++) {
            buffer[latentIndices[i]] = latentDraw[i];
        }

        return buffer;
    }


}

//    public static void safeSolveSymmPosDef(DenseMatrix64F A,
//                                           WrappedVector b,
//                                           WrappedVector x) {
//        final int dim = b.getDim();
//
//        assert (A.getNumRows() == dim && A.getNumCols() == dim);
//
//        final DenseMatrix64F B = wrap(b.getBuffer(), b.getOffset(), dim, 1);
//        final DenseMatrix64F X = new DenseMatrix64F(dim, 1);
//
//        safeSolveSymmPosDef(A, B, X);
//
//
//        for (int row = 0; row < dim; ++row) {
//            x.set(row, X.unsafe_get(row, 0));
//        }
//    }
//
//    public static void safeSolveSymmPosDef(DenseMatrix64F A, DenseMatrix64F B, DenseMatrix64F X) {
//
//        final int finiteCount = countFiniteNonZeroDiagonals(A);
//
//        InversionResult result;
//        if (finiteCount == 0) {
//            Arrays.fill(X.getData(), 0);
//        } else {
//            LinearSolver<DenseMatrix64F> solver = LinearSolverFactory.symmPosDef(A.getNumCols());
//            DenseMatrix64F Abis = new DenseMatrix64F(A);
//            if(solver.setA(Abis)) {
//                solver.solve(B, X);
//            } else {
//                LinearSolver<DenseMatrix64F> solverSVD = LinearSolverFactory.pseudoInverse(true);
//                solverSVD.setA(A);
//                solverSVD.solve(B, X);
//            }
//        }
//    }


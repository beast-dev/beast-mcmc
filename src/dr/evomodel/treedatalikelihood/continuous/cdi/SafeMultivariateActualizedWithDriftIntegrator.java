package dr.evomodel.treedatalikelihood.continuous.cdi;

import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

import static dr.math.matrixAlgebra.missingData.MissingOps.*;

/**
 * @author Marc A. Suchard
 * @author Paul Bastide
 */

public class SafeMultivariateActualizedWithDriftIntegrator extends SafeMultivariateDiagonalActualizedWithDriftIntegrator {

    private static boolean DEBUG = false;

    public SafeMultivariateActualizedWithDriftIntegrator(PrecisionType precisionType, int numTraits, int dimTrait, int bufferCount,
                                                         int diffusionCount) {
        super(precisionType, numTraits, dimTrait, bufferCount, diffusionCount);

        allocateStorage();

        System.err.println("Trying SafeMultivariateActualizedWithDriftIntegrator");
    }

    @Override
    public void getBranchActualization(int bufferIndex, double[] actualization) {

        if (bufferIndex == -1) {
            throw new RuntimeException("Not yet implemented");
        }

        assert (actualization != null);
        assert (actualization.length >= dimTrait * dimTrait);

        System.arraycopy(actualizations, bufferIndex * dimTrait * dimTrait,
                actualization, 0, dimTrait * dimTrait);
    }

    private static final boolean TIMING = false;

    private void allocateStorage() {

        actualizations = new double[dimTrait * dimTrait * bufferCount];

        matrixQdiPip = new DenseMatrix64F(dimTrait, dimTrait);
        matrixQdjPjp = new DenseMatrix64F(dimTrait, dimTrait);
    }

    @Override
    public void setDiffusionStationaryVariance(int precisionIndex, final double[] alphaEig, final double[] alphaRot) {

        super.setDiffusionStationaryVariance(precisionIndex, alphaEig, alphaRot);

        // Transform back in original space
        final int offset = dimTrait * dimTrait * precisionIndex;
        DenseMatrix64F stationaryVariance = wrap(stationaryVariances, offset, dimTrait, dimTrait);
//        DenseMatrix64F rot = wrap(alphaRot, 0, dimTrait, dimTrait);
        transformMatrixBack(stationaryVariances, offset, alphaRot, 0);

        if (DEBUG) {
            System.err.println("At precision index: " + precisionIndex);
//            System.err.println("variance : " + variance);
            System.err.println("stationary variance: " + stationaryVariance);
        }
    }

    private void transformMatrixBack(double[] matrixDouble, int matrixOffset, double[] rotationDouble, int rotationOffset) {
        DenseMatrix64F matrix = wrap(matrixDouble, matrixOffset, dimTrait, dimTrait);
        DenseMatrix64F rotation = wrap(rotationDouble, rotationOffset, dimTrait, dimTrait);
        transformMatrixBack(matrix, rotation);
        unwrap(matrix, matrixDouble, matrixOffset);
    }

    private void transformDiagonalMatrixBack(double[] diagonalMatrix, double[] matrixDestination, int matrixOffset, double[] rotationDouble, int rotationOffset) {
        DenseMatrix64F matrix = wrapDiagonal(diagonalMatrix, matrixOffset, dimTrait);
        DenseMatrix64F rotation = wrap(rotationDouble, rotationOffset, dimTrait, dimTrait);
        transformMatrixBack(matrix, rotation);
        unwrap(matrix, matrixDestination, matrixOffset);
    }

    private void transformMatrixBack(DenseMatrix64F matrix, DenseMatrix64F rotation) {
        DenseMatrix64F tmp = new DenseMatrix64F(dimTrait, dimTrait);
        CommonOps.multTransB(matrix, rotation, tmp);
        CommonOps.mult(rotation, tmp, matrix);
    }

//    private void computeStationaryVariance(DenseMatrix64F diffusion, DenseMatrix64F alphaMatrix, DenseMatrix64F stationaryVariance){
//        DenseMatrix64F variance_vec = vectorize(diffusion);
//        DenseMatrix64F kro_sum_A_inv = sumKronecker(alphaMatrix, alphaMatrix);
//        CommonOps.invert(kro_sum_A_inv);
//        DenseMatrix64F stationaryVarianceVec = new DenseMatrix64F(dimTrait * dimTrait, 1);
//        CommonOps.mult(kro_sum_A_inv, variance_vec, stationaryVarianceVec);
//        unVectorizeSquare(stationaryVarianceVec, stationaryVariance);
//    }

//    private DenseMatrix64F vectorize(DenseMatrix64F A){
//        int m = A.numRows;
//        int n = A.numCols;
//        DenseMatrix64F B = new DenseMatrix64F(m * n, 1);
//
//        for (int i = 0; i < n; ++i){
//            CommonOps.extract(A, 0, m, i, i+1, B, i * m, 0);
//        }
//
//        return B;
//    }

//    private DenseMatrix64F sumKronecker(DenseMatrix64F A, DenseMatrix64F B){
//        int m = A.numCols;
//        int n = B.numCols;
//
//        if (m != A.numRows || n != B.numRows){
//            throw new RuntimeException("Wrong dimensions in Kronecker sum");
//        }
//
//        DenseMatrix64F C1 = new DenseMatrix64F(m * n, m * n);
//        DenseMatrix64F C2 = new DenseMatrix64F(m * n, m * n);
//
//        DenseMatrix64F I_m = CommonOps.identity(m);
//        DenseMatrix64F I_n = CommonOps.identity(n);
//
//        CommonOps.kron(A, I_n, C1);
//        CommonOps.kron(I_m, B, C2);
//        CommonOps.addEquals(C1, C2);
//
//        return C1;
//    }

//    private void unVectorizeSquare(DenseMatrix64F vector, DenseMatrix64F matrix){
//        int n = matrix.numRows;
//        if (1 != vector.numCols || n * n != vector.numRows || n != matrix.numCols){
//            throw new RuntimeException("Wrong dimensions in unVectorizeSquare");
//        }
//
//        for (int i = 0; i < n; ++i){
//            CommonOps.extract(vector, i * n, (i+1) * n, 0, 1, matrix, 0, i);
//        }
//    }

    @Override
    public void updateOrnsteinUhlenbeckDiffusionMatrices(int precisionIndex, final int[] probabilityIndices,
                                                         final double[] edgeLengths, final double[] optimalRates,
                                                         final double[] diagonalStrengthOfSelectionMatrix,
                                                         final double[] rotation,
                                                         int updateCount) {

        assert (diffusions != null);
        assert (probabilityIndices.length >= updateCount);
        assert (edgeLengths.length >= updateCount);

        if (DEBUG) {
            System.err.println("Matrices (safe with drift):");
        }

        final int matrixSize = dimTrait * dimTrait;
        final int unscaledOffset = matrixSize * precisionIndex;


        if (TIMING) {
            startTime("actualization");
        }

        for (int up = 0; up < updateCount; ++up) {

            if (DEBUG) {
                System.err.println("\t" + probabilityIndices[up] + " <- " + edgeLengths[up]);
            }

            final double edgeLength = edgeLengths[up];
            branchLengths[dimMatrix * probabilityIndices[up]] = edgeLength;  // TODO Remove dimMatrix

            final int scaledOffset = matrixSize * probabilityIndices[up];


            double[] diagonalActualizations = new double[matrixSize];
            computeDiagonalActualization(diagonalStrengthOfSelectionMatrix, edgeLength, dimTrait, diagonalActualizations, 0);
            transformDiagonalMatrixBack(diagonalActualizations, actualizations, scaledOffset, rotation, 0);
        }

        if (TIMING) {
            endTime("actualization");
        }


        if (TIMING) {
            startTime("diffusion");
        }

        for (int up = 0; up < updateCount; ++up) {

            if (DEBUG) {
                System.err.println("\t" + probabilityIndices[up] + " <- " + edgeLengths[up]);
            }

            final int scaledOffset = matrixSize * probabilityIndices[up];

            computeVarianceBranch(stationaryVariances, unscaledOffset, dimTrait, actualizations, variances, scaledOffset);
            invertVector(variances, precisions, scaledOffset, dimTrait);
        }

        if (TIMING) {
            endTime("diffusion");
        }

        assert (optimalRates != null);
        assert (displacements != null);
        assert (optimalRates.length >= updateCount * dimTrait);

        if (TIMING) {
            startTime("drift1");
        }

        int offset = 0;
        for (int up = 0; up < updateCount; ++up) {

            final int scaledOffset = matrixSize * probabilityIndices[up];
            final int pio = dimTrait * probabilityIndices[up];

            computeDisplacement(optimalRates, offset, actualizations, scaledOffset, displacements, pio, dimTrait);
            offset += dimTrait;
        }

        if (TIMING) {
            endTime("drift1");
        }

        precisionOffset = dimTrait * dimTrait * precisionIndex;
        precisionLogDet = determinants[precisionIndex];
    }

//    private static void computeActualization(final double[] source,
//                                             final double edgeLength,
//                                             final int dim,
//                                             final double[] destination,
//                                             final int destinationOffset) {
//        DenseMatrix64F alphaMatrix = wrap(source, 0, dim, dim);
//        DenseMatrix64F actualization = new DenseMatrix64F(dim, dim);
//        scaledMatrixExponential(alphaMatrix, -edgeLength, actualization); // QUESTION: Does this already exist ?
//        unwrap(actualization, destination, destinationOffset);
//    }

    private static void computeVarianceBranch(final double[] source,
                                              final int sourceOffset,
                                              final int dim,
                                              final double[] actualizations,
                                              final double[] destination,
                                              final int destinationOffset) {
        DenseMatrix64F actualization = wrap(actualizations, destinationOffset, dim, dim);
        DenseMatrix64F variance = wrap(source, sourceOffset, dim, dim);
        DenseMatrix64F temp = new DenseMatrix64F(dim, dim);

        CommonOps.multTransA(variance, actualization, temp);
        CommonOps.multAdd(-1.0, actualization, temp, variance);

        unwrap(variance, destination, destinationOffset);
    }

    private static void invertVector(final double[] source,
                                     final double[] destination,
                                     final int offset,
                                     final int dim) {
        DenseMatrix64F sourceMatrix = wrap(source, offset, dim, dim);
        DenseMatrix64F destinationMatrix = new DenseMatrix64F(dim, dim);

        CommonOps.invert(sourceMatrix, destinationMatrix);

        unwrap(destinationMatrix, destination, offset);
    }


    private static void computeDisplacement(final double[] source,
                                            final int sourceOffset,
                                            final double[] actualizations,
                                            final int actualizationOffset,
                                            final double[] destination,
                                            final int destinationOffset,
                                            final int dim) {
        DenseMatrix64F actualization = wrap(actualizations, actualizationOffset, dim, dim);
        DenseMatrix64F optVal = wrap(source, sourceOffset, dim, 1);
        DenseMatrix64F temp = CommonOps.identity(dim);
        DenseMatrix64F displacement = new DenseMatrix64F(dim, 1);

        CommonOps.addEquals(temp, -1.0, actualization);
        CommonOps.mult(temp, optVal, displacement);

        unwrap(displacement, destination, destinationOffset);
    }

//    private static void eigenStrengthOfSelectionMatrix(DenseMatrix64F A, EigenDecomposition eigA){
//        int n = A.numCols;
//        if (n != A.numRows) throw new RuntimeException("Selection strength A matrix must be square.");
////        EigenDecomposition eigA = DecompositionFactory.eig(n, true);
//        if( !eigA.decompose(A) ) throw new RuntimeException("Eigen decomposition failed.");
//        for (int p = 0; p < n; ++p) {
//            if (!eigA.getEigenvalue(p).isReal()) throw new RuntimeException("Selection strength A should only have real eigenvalues.");
//        }
//    }

//    private static void scaledMatrixExponential(DenseMatrix64F A, double lambda, DenseMatrix64F C){
//        int n = A.numCols;
//        if (n != A.numRows) throw new RuntimeException("Selection strength A matrix must be square.");
//        EigenDecomposition eigA = DecompositionFactory.eig(n, true);
//        if( !eigA.decompose(A) ) throw new RuntimeException("Eigen decomposition failed.");
//        DenseMatrix64F expDiag = CommonOps.identity(n);
//        for (int p = 0; p < n; ++p) {
//            Complex64F ev = eigA.getEigenvalue(p);
//            if (!ev.isReal()) throw new RuntimeException("Selection strength A should only have real eigenvalues.");
//            expDiag.set(p, p, Math.exp(lambda * ev.real));
//        }
//        DenseMatrix64F V = EigenOps.createMatrixV(eigA);
//        DenseMatrix64F tmp = new DenseMatrix64F(n, n);
//        CommonOps.mult(V, expDiag, tmp);
//        CommonOps.invert(V);
//        CommonOps.mult(tmp, V, C);
//    }

    @Override
    public void updatePreOrderPartial(
            final int kBuffer, // parent
            final int iBuffer, // node
            final int iMatrix,
            final int jBuffer, // sibling
            final int jMatrix) {

        throw new RuntimeException("Not yet implemented");
    }

    private double[] actualizations;
    private DenseMatrix64F matrixQdiPip;
    private DenseMatrix64F matrixQdjPjp;

    @Override
    void computePartialPrecision(int ido, int jdo, int imo, int jmo,
                                 DenseMatrix64F Pip, DenseMatrix64F Pjp, DenseMatrix64F Pk) {

        final DenseMatrix64F Qdi = wrap(actualizations, imo, dimTrait, dimTrait);
        final DenseMatrix64F Qdj = wrap(actualizations, jmo, dimTrait, dimTrait);

        final DenseMatrix64F QdiPip = matrixQdiPip;
        CommonOps.multTransA(Qdi, Pip, QdiPip);

        final DenseMatrix64F QdjPjp = matrixQdjPjp;
        CommonOps.multTransA(Qdj, Pjp, QdjPjp);

        final DenseMatrix64F QdiPipQdi = matrix0;
        final DenseMatrix64F QdjPjpQdj = matrix1;
        CommonOps.mult(QdiPip, Qdi, QdiPipQdi);
        CommonOps.mult(QdjPjp, Qdj, QdjPjpQdj);
        CommonOps.add(QdiPipQdi, QdjPjpQdj, Pk);

        if (DEBUG) {
            System.err.println("Qdi: " + Qdi);
            System.err.println("\tQdiPip: " + QdiPip);
            System.err.println("\tQdiPipQdi: " + QdiPipQdi);
            System.err.println("\tQdj: " + Qdj);
            System.err.println("\tQdjPjp: " + QdjPjp);
            System.err.println("\tQdjPjpQdj: " + QdjPjpQdj);
        }
    }

    @Override
    void computeWeightedSum(final double[] ipartial,
                            final double[] jpartial,
                            final int dimTrait,
                            final double[] out) {
        weightedSum(ipartial, 0, matrixQdiPip, jpartial, 0, matrixQdjPjp, dimTrait, out);
    }

//    public double[] getStationaryVariance(int precisionIndex) {
//
//        assert (stationaryVariances != null);
//
//        final int offset = dimTrait * dimTrait * precisionIndex;
//
//        double[] buffer = new double[dimTrait * dimTrait];
//
//        System.arraycopy(stationaryVariances, offset, buffer, 0, dimTrait * dimTrait);
//
//        return buffer;
//    }
}

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

    public SafeMultivariateActualizedWithDriftIntegrator(PrecisionType precisionType,
                                                         int numTraits, int dimTrait, int bufferCount,
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

//    private static final boolean TIMING = false;

    private void allocateStorage() {

        actualizations = new double[dimTrait * dimTrait * bufferCount];

        matrixQdiPip = new DenseMatrix64F(dimTrait, dimTrait);
        matrixQdjPjp = new DenseMatrix64F(dimTrait, dimTrait);
    }

    ///////////////////////////////////////////////////////////////////////////
    /// Setting variances, displacement and actualization vectors
    ///////////////////////////////////////////////////////////////////////////

    @Override
    public void setDiffusionStationaryVariance(int precisionIndex, final double[] alphaEig, final double[] alphaRot) {

        super.setDiffusionStationaryVariance(precisionIndex, alphaEig, alphaRot);

        // Transform back in original space
        final int offset = dimTrait * dimTrait * precisionIndex;
        DenseMatrix64F stationaryVariance = wrap(stationaryVariances, offset, dimTrait, dimTrait);
        transformMatrixBack(stationaryVariances, offset, alphaRot, 0);

        if (DEBUG) {
            System.err.println("At precision index: " + precisionIndex);
            System.err.println("stationary variance: " + stationaryVariance);
        }
    }

    @Override
    void setStationaryVariance(int offset, double[] scales, int matrixSize, double[] rotation) {
        assert (rotation.length == matrixSize);
        DenseMatrix64F rotMat = wrap(rotation, 0, dimTrait, dimTrait);
        DenseMatrix64F variance = wrap(inverseDiffusions, offset, dimTrait, dimTrait);
        transformMatrix(variance, rotMat);
        double[] transVar = new double[matrixSize];
        unwrap(variance, transVar, 0);
        scaleInv(transVar, 0, scales, stationaryVariances, offset, matrixSize);
    }

//    private void transformMatrix(DenseMatrix64F matrix, DenseMatrix64F rotation) {
//        DenseMatrix64F tmp = new DenseMatrix64F(dimTrait, dimTrait);
//        CommonOps.invert(rotation); // Warning: side effect on rotation matrix.
//        CommonOps.mult(rotation, matrix, tmp);
//        CommonOps.multTransB(tmp, rotation, matrix);
//    }

    private void transformMatrix(DenseMatrix64F matrix, DenseMatrix64F rotation) {
        DenseMatrix64F tmp = new DenseMatrix64F(dimTrait, dimTrait);
        CommonOps.multTransA(rotation, matrix, tmp); // TODO: this is a specialized version for symmetric A (see above general version)
        CommonOps.mult(tmp, rotation, matrix);
    }

    private void transformMatrixBack(double[] matrixDouble, int matrixOffset,
                                     double[] rotationDouble, int rotationOffset) {
        DenseMatrix64F matrix = wrap(matrixDouble, matrixOffset, dimTrait, dimTrait);
        DenseMatrix64F rotation = wrap(rotationDouble, rotationOffset, dimTrait, dimTrait);
        transformMatrixBack(matrix, rotation);
        unwrap(matrix, matrixDouble, matrixOffset);
    }

    private void transformDiagonalMatrixBack(double[] diagonalMatrix, double[] matrixDestination, int matrixOffset,
                                             double[] rotationDouble, int rotationOffset) {
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

    @Override
    void computeActualization(final double[] diagonalStrengthOfSelectionMatrix,
                              final double[] rotation,
                              final double edgeLength,
                              final int scaledOffsetDiagonal,
                              final int scaledOffset) {
        double[] diagonalActualizations = new double[dimTrait * dimTrait];
        computeDiagonalActualization(diagonalStrengthOfSelectionMatrix, edgeLength, dimTrait,
                diagonalActualizations, 0);
        transformDiagonalMatrixBack(diagonalActualizations, actualizations, scaledOffset, rotation, 0);
    }

    @Override
    void computeVarianceBranch(final int sourceOffset,
                               final int destinationOffset,
                               final int destinationOffsetDiagonal) {
        DenseMatrix64F actualization = wrap(actualizations, destinationOffset, dimTrait, dimTrait);
        DenseMatrix64F variance = wrap(stationaryVariances, sourceOffset, dimTrait, dimTrait);
        DenseMatrix64F temp = new DenseMatrix64F(dimTrait, dimTrait);

        CommonOps.multTransA(variance, actualization, temp);
        CommonOps.multAdd(-1.0, actualization, temp, variance);

        unwrap(variance, variances, destinationOffset);
    }

    @Override
    void computeActualizedDisplacement(final double[] optimalRates,
                                       final int offset,
                                       final int up,
                                       final int actualizationOffset,
                                       final int pio) {
        DenseMatrix64F actualization = wrap(actualizations, actualizationOffset, dimTrait, dimTrait);
        DenseMatrix64F optVal = wrap(optimalRates, offset, dimTrait, 1);
        DenseMatrix64F temp = CommonOps.identity(dimTrait);
        DenseMatrix64F displacement = new DenseMatrix64F(dimTrait, 1);

        CommonOps.addEquals(temp, -1.0, actualization);
        CommonOps.mult(temp, optVal, displacement);

        unwrap(displacement, displacements, pio);
    }

    ///////////////////////////////////////////////////////////////////////////
    /// Tree-traversal functions
    ///////////////////////////////////////////////////////////////////////////

    @Override
    public void updatePreOrderPartial(
            final int kBuffer, // parent
            final int iBuffer, // node
            final int iMatrix,
            final int jBuffer, // sibling
            final int jMatrix) {

        throw new RuntimeException("Not yet implemented");
    }

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

    private double[] actualizations;
    private DenseMatrix64F matrixQdiPip;
    private DenseMatrix64F matrixQdjPjp;
}

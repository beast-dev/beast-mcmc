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
                                                         int diffusionCount,
                                                         boolean isActualizationSymmetric) {
        super(precisionType, numTraits, dimTrait, bufferCount, diffusionCount);

        allocateStorage();

        this.isActualizationSymmetric = isActualizationSymmetric;

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

    @Override
    public void getBranchExpectation(double[] actualization, double[] parentValue, double[] displacement,
                                     double[] expectation) {

        assert (expectation != null);
        assert (expectation.length >= dimTrait);

        assert (actualization != null);
        assert (actualization.length >= dimTrait * dimTrait);

        assert (parentValue != null);
        assert (parentValue.length >= dimTrait);

        assert (displacement != null);
        assert (displacement.length >= dimTrait);

        DenseMatrix64F branchExpectationMatrix = new DenseMatrix64F(dimTrait, 1);
        CommonOps.mult(wrap(actualization, 0, dimTrait, dimTrait),
                wrap(parentValue, 0, dimTrait, 1),
                branchExpectationMatrix);
        CommonOps.addEquals(branchExpectationMatrix, wrap(displacement, 0, dimTrait, 1));

        unwrap(branchExpectationMatrix, expectation, 0);
    }

//    private static final boolean TIMING = false;

    private void allocateStorage() {

        actualizations = new double[dimTrait * dimTrait * bufferCount];

        matrixQdiPip = new DenseMatrix64F(dimTrait, dimTrait);
        matrixQdjPjp = new DenseMatrix64F(dimTrait, dimTrait);

        matrixNiacc = new DenseMatrix64F(dimTrait, 1);
    }

    ///////////////////////////////////////////////////////////////////////////
    /// Setting variances, displacement and actualization vectors
    ///////////////////////////////////////////////////////////////////////////

    @Override
    public void setDiffusionStationaryVariance(int precisionIndex, final double[] alphaEig, final double[] alphaRot) {

        super.setDiffusionStationaryVariance(precisionIndex, alphaEig, alphaRot);

        // Transform back in original space
        final int offset = dimTrait * dimTrait * precisionIndex;
        transformMatrixBack(stationaryVariances, offset, alphaRot, 0);

        if (DEBUG) {
            System.err.println("At precision index: " + precisionIndex);
            System.err.println("stationary variance: " + wrap(stationaryVariances, offset, dimTrait, dimTrait));
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

    private void transformMatrix(DenseMatrix64F matrix, DenseMatrix64F rotation) {
        if (isActualizationSymmetric) {
            transformMatrixSymmetric(matrix, rotation);
        } else {
            transformMatrixGeneral(matrix, rotation);
        }
    }

    private void transformMatrixGeneral(DenseMatrix64F matrix, DenseMatrix64F rotation) {
        DenseMatrix64F tmp = new DenseMatrix64F(dimTrait, dimTrait);
        CommonOps.invert(rotation); // Warning: side effect on rotation matrix.
        CommonOps.mult(rotation, matrix, tmp);
        CommonOps.multTransB(tmp, rotation, matrix);
    }

    private void transformMatrixSymmetric(DenseMatrix64F matrix, DenseMatrix64F rotation) {
        DenseMatrix64F tmp = new DenseMatrix64F(dimTrait, dimTrait);
        CommonOps.multTransA(rotation, matrix, tmp);
        CommonOps.mult(tmp, rotation, matrix);
    }

    private void transformMatrixBack(double[] matrixDouble, int matrixOffset,
                                     double[] rotationDouble, int rotationOffset) {
        DenseMatrix64F matrix = wrap(matrixDouble, matrixOffset, dimTrait, dimTrait);
        DenseMatrix64F rotation = wrap(rotationDouble, rotationOffset, dimTrait, dimTrait);
        transformMatrixBack(matrix, rotation);
        unwrap(matrix, matrixDouble, matrixOffset);
    }

    private void transformMatrixBack(DenseMatrix64F matrix, DenseMatrix64F rotation) {
        DenseMatrix64F tmp = new DenseMatrix64F(dimTrait, dimTrait);
        CommonOps.multTransB(matrix, rotation, tmp);
        CommonOps.mult(rotation, tmp, matrix);
    }

    private void transformDiagonalMatrixBack(double[] diagonalMatrix, double[] matrixDestination, int matrixOffset,
                                             double[] rotationDouble, int rotationOffset) {
        DenseMatrix64F matrix = wrapDiagonal(diagonalMatrix, matrixOffset, dimTrait);
        DenseMatrix64F rotation = wrap(rotationDouble, rotationOffset, dimTrait, dimTrait);
        transformMatrixBase(matrix, rotation);
        unwrap(matrix, matrixDestination, matrixOffset);
    }

    private void transformMatrixBase(DenseMatrix64F matrix, DenseMatrix64F rotation) {
        if (isActualizationSymmetric) {
            transformMatrixBack(matrix, rotation);
        } else {
            transformMatrixBaseGeneral(matrix, rotation);
        }
    }

    private void transformMatrixBaseGeneral(DenseMatrix64F matrix, DenseMatrix64F rotation) {
        DenseMatrix64F tmp = new DenseMatrix64F(dimTrait, dimTrait);
        CommonOps.mult(rotation, matrix, tmp);
        CommonOps.invert(rotation); // Warning: side effect on rotation matrix.
        CommonOps.mult(tmp, rotation, matrix);
    }

    @Override
    void computeActualization(final double[] diagonalStrengthOfSelectionMatrix,
                              final double[] rotation,
                              final double edgeLength,
                              final int scaledOffsetDiagonal,
                              final int scaledOffset) {
        if (edgeLength == 0.0) {
            unwrapIdentity(actualizations, scaledOffset, dimTrait);
        } else {
            double[] diagonalActualizations = new double[dimTrait * dimTrait];
            computeDiagonalActualization(diagonalStrengthOfSelectionMatrix, edgeLength, dimTrait,
                    diagonalActualizations, 0);
            transformDiagonalMatrixBack(diagonalActualizations, actualizations, scaledOffset, rotation, 0);
        }
    }

    @Override
    void computeVarianceBranch(final int sourceOffset,
                               final int destinationOffset,
                               final int destinationOffsetDiagonal,
                               final double edgeLength) {
        DenseMatrix64F actualization = wrap(actualizations, destinationOffset, dimTrait, dimTrait);
        DenseMatrix64F variance = wrap(stationaryVariances, sourceOffset, dimTrait, dimTrait);
        DenseMatrix64F temp = new DenseMatrix64F(dimTrait, dimTrait);

        CommonOps.multTransB(variance, actualization, temp);
        CommonOps.multAdd(-1.0, actualization, temp, variance);

        unwrap(variance, variances, destinationOffset);
    }

    @Override
    void computeActualizedDisplacement(final double[] optimalRates,
                                       final int offset,
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

//    @Override
//    public void updatePreOrderPartial(
//            final int kBuffer, // parent
//            final int iBuffer, // node
//            final int iMatrix,
//            final int jBuffer, // sibling
//            final int jMatrix) {
//
//        throw new RuntimeException("Not yet implemented");
//    }

    @Override
    void actualizePrecision(DenseMatrix64F Pjp, DenseMatrix64F QjPjp, int jbo, int jmo, int jdo) {
        final DenseMatrix64F Qdj = wrap(actualizations, jmo, dimTrait, dimTrait);
        scalePrecision(Qdj, Pjp, QjPjp, Pjp);
    }

    @Override
    void actualizeVariance(DenseMatrix64F Vip, int ibo, int imo, int ido) {
        final DenseMatrix64F Qdi = wrap(actualizations, imo, dimTrait, dimTrait);
        final DenseMatrix64F QiVip = matrixQdiPip;
        scaleVariance(Qdi, Vip, QiVip, Vip);
    }

    @Override
    void scaleAndDriftMean(int ibo, int imo, int ido) {
        final DenseMatrix64F Qdi = wrap(actualizations, imo, dimTrait, dimTrait);
        final DenseMatrix64F ni = wrap(preOrderPartials, ibo, dimTrait, 1);
        final DenseMatrix64F niacc = matrixNiacc;
        CommonOps.mult(Qdi, ni, niacc);
        unwrap(niacc, preOrderPartials, ibo);

        for (int g = 0; g < dimTrait; ++g) {
            preOrderPartials[ibo + g] += displacements[ido + g];
        }

    }

    @Override
    void computePartialPrecision(int ido, int jdo, int imo, int jmo,
                                 DenseMatrix64F Pip, DenseMatrix64F Pjp, DenseMatrix64F Pk) {

        final DenseMatrix64F Qdi = wrap(actualizations, imo, dimTrait, dimTrait);
        final DenseMatrix64F Qdj = wrap(actualizations, jmo, dimTrait, dimTrait);

        final DenseMatrix64F QdiPip = matrixQdiPip;
        final DenseMatrix64F QdiPipQdi = matrix0;
        scalePrecision(Qdi, Pip, QdiPip, QdiPipQdi);

        final DenseMatrix64F QdjPjpQdj = matrix1;
        final DenseMatrix64F QdjPjp = matrixQdjPjp;
        scalePrecision(Qdj, Pjp, QdjPjp, QdjPjpQdj);

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

    private void scalePrecision(DenseMatrix64F Q, DenseMatrix64F P,
                                DenseMatrix64F QtP, DenseMatrix64F QtPQ) {
        CommonOps.multTransA(Q, P, QtP);
        CommonOps.mult(QtP, Q, QtPQ);
    }

    private void scaleVariance(DenseMatrix64F Q, DenseMatrix64F P,
                               DenseMatrix64F QtP, DenseMatrix64F QtPQ) {
        CommonOps.mult(Q, P, QtP);
        CommonOps.multTransB(QtP, Q, QtPQ);
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
    private DenseMatrix64F matrixNiacc;
    private final boolean isActualizationSymmetric;
}

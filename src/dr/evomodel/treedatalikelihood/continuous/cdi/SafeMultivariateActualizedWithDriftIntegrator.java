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
                                                         int numTraits, int dimTrait, int dimProcess,
                                                         int bufferCount, int diffusionCount,
                                                         boolean isActualizationSymmetric) {
        super(precisionType, numTraits, dimTrait, dimProcess, bufferCount, diffusionCount);

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

    private static final boolean TIMING = false;

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

        int dim = alphaEig.length;
        assert alphaRot.length == dim * dim;

        super.setDiffusionStationaryVariance(precisionIndex, alphaEig, alphaRot);

        // Transform back in original space
        final int offset = dimProcess * dimProcess * precisionIndex;
        transformMatrixBack(stationaryVariances, offset, alphaRot, 0);

        if (DEBUG) {
            System.err.println("At precision index: " + precisionIndex);
            System.err.println("stationary variance: " + wrap(stationaryVariances, offset, dim, dim));
        }
    }

    @Override
    void setStationaryVariance(int offset, double[] scales, int matrixSize, double[] rotation) {
        assert (rotation.length == matrixSize);

        DenseMatrix64F rotMat = wrap(rotation, 0, dimProcess, dimProcess);
        DenseMatrix64F variance = wrap(inverseDiffusions, offset, dimProcess, dimProcess);
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
        DenseMatrix64F tmp = new DenseMatrix64F(dimProcess, dimProcess);
        CommonOps.invert(rotation); // Warning: side effect on rotation matrix.
        CommonOps.mult(rotation, matrix, tmp);
        CommonOps.multTransB(tmp, rotation, matrix);
    }

    private void transformMatrixSymmetric(DenseMatrix64F matrix, DenseMatrix64F rotation) {
        DenseMatrix64F tmp = new DenseMatrix64F(dimProcess, dimProcess);
        CommonOps.multTransA(rotation, matrix, tmp);
        CommonOps.mult(tmp, rotation, matrix);
    }

    private void transformMatrixBack(double[] matrixDouble, int matrixOffset,
                                     double[] rotationDouble, int rotationOffset) {
        DenseMatrix64F matrix = wrap(matrixDouble, matrixOffset, dimProcess, dimProcess);
        DenseMatrix64F rotation = wrap(rotationDouble, rotationOffset, dimProcess, dimProcess);
        transformMatrixBack(matrix, rotation);
        unwrap(matrix, matrixDouble, matrixOffset);
    }

    private void transformMatrixBack(DenseMatrix64F matrix, DenseMatrix64F rotation) {
        DenseMatrix64F tmp = new DenseMatrix64F(dimProcess, dimProcess);
        CommonOps.multTransB(matrix, rotation, tmp);
        CommonOps.mult(rotation, tmp, matrix);
    }

    private void transformDiagonalMatrixBack(double[] diagonalMatrix, double[] matrixDestination, int matrixOffset,
                                             double[] rotationDouble, int rotationOffset) {
        DenseMatrix64F matrix = wrapDiagonal(diagonalMatrix, matrixOffset, dimProcess);
        DenseMatrix64F rotation = wrap(rotationDouble, rotationOffset, dimProcess, dimProcess);
        transformMatrixBase(matrix, rotation);
        unwrap(matrix, matrixDestination, matrixOffset);
    }

    private DenseMatrix64F getInverseSelectionStrength(double[] diagonalMatrix, double[] rotationDouble) {
        DenseMatrix64F matrix = wrapDiagonalInverse(diagonalMatrix, 0, dimProcess);
        DenseMatrix64F rotation = wrap(rotationDouble, 0, dimProcess, dimProcess);
        transformMatrixBase(matrix, rotation);
        return matrix;
    }

    private void transformMatrixBase(DenseMatrix64F matrix, DenseMatrix64F rotation) {
        if (isActualizationSymmetric) {
            transformMatrixBack(matrix, rotation);
        } else {
            transformMatrixBaseGeneral(matrix, rotation);
        }
    }

    private void transformMatrixBaseGeneral(DenseMatrix64F matrix, DenseMatrix64F rotation) {
        DenseMatrix64F tmp = new DenseMatrix64F(dimProcess, dimProcess);
        CommonOps.mult(rotation, matrix, tmp);
        CommonOps.invert(rotation); // Warning: side effect on rotation matrix.
        CommonOps.mult(tmp, rotation, matrix);
    }

    @Override
    void computeOUActualization(final double[] diagonalStrengthOfSelectionMatrix,
                                final double[] rotation,
                                final double edgeLength,
                                final int scaledOffsetDiagonal,
                                final int scaledOffset) {
        double[] diagonalActualizations = new double[dimTrait];
        computeOUDiagonalActualization(diagonalStrengthOfSelectionMatrix, edgeLength, dimProcess,
                diagonalActualizations, 0);
        transformDiagonalMatrixBack(diagonalActualizations, actualizations, scaledOffset, rotation, 0);
    }

    @Override
    void computeOUVarianceBranch(final int sourceOffset,
                                 final int destinationOffset,
                                 final int destinationOffsetDiagonal,
                                 final double edgeLength) {
        DenseMatrix64F actualization = wrap(actualizations, destinationOffset, dimProcess, dimProcess);
        DenseMatrix64F variance = wrap(stationaryVariances, sourceOffset, dimProcess, dimProcess);
        DenseMatrix64F temp = new DenseMatrix64F(dimProcess, dimProcess);

        CommonOps.multTransB(variance, actualization, temp);
        CommonOps.multAdd(-1.0, actualization, temp, variance);

        unwrap(variance, variances, destinationOffset);
    }

    @Override
    void computeOUActualizedDisplacement(final double[] optimalRates,
                                         final int offset,
                                         final int actualizationOffset,
                                         final int pio) {
        DenseMatrix64F actualization = wrap(actualizations, actualizationOffset, dimProcess, dimProcess);
        DenseMatrix64F optVal = wrap(optimalRates, offset, dimProcess, 1);
        DenseMatrix64F temp = CommonOps.identity(dimProcess);
        DenseMatrix64F displacement = new DenseMatrix64F(dimProcess, 1);

        CommonOps.addEquals(temp, -1.0, actualization);
        CommonOps.mult(temp, optVal, displacement);

        unwrap(displacement, displacements, pio);
    }

    private void computeIOUActualizedDisplacement(final double[] optimalRates,
                                                  final int offset,
                                                  final int pio,
                                                  double branchLength,
                                                  DenseMatrix64F inverseSelectionStrength) {
        DenseMatrix64F displacementOU = wrap(displacements, pio, dimProcess, 1);
        DenseMatrix64F optVal = wrap(optimalRates, offset, dimProcess, 1);
        DenseMatrix64F displacement = new DenseMatrix64F(dimProcess, 1);

        CommonOps.mult(inverseSelectionStrength, displacementOU, displacement);
        CommonOps.scale(-1.0, displacement);

        CommonOps.addEquals(displacement, branchLength, optVal);

        unwrap(displacement, displacements, pio + dimProcess);
    }

    private void computeIOUVarianceBranch(final int sourceOffset,
                                          final int destinationOffset,
                                          double branchLength,
                                          DenseMatrix64F inverseSelectionStrength) {
        DenseMatrix64F actualization = wrap(actualizations, destinationOffset, dimProcess, dimProcess);
        DenseMatrix64F stationaryVariance = wrap(stationaryVariances, sourceOffset, dimProcess, dimProcess);

        DenseMatrix64F invAS = new DenseMatrix64F(dimProcess, dimProcess);
        CommonOps.mult(inverseSelectionStrength, stationaryVariance, invAS);

        //// Variance YY
        DenseMatrix64F varianceYY = wrap(variances, destinationOffset, dimProcess, dimProcess);

        //// Variance XX
        DenseMatrix64F varianceXX = new DenseMatrix64F(dimProcess, dimProcess);
        // Variance 1
        CommonOps.multTransB(invAS, inverseSelectionStrength, varianceXX);
        DenseMatrix64F temp = new DenseMatrix64F(dimProcess, dimProcess);
        CommonOps.multTransB(varianceXX, actualization, temp);
        CommonOps.multAdd(-1.0, actualization, temp, varianceXX);
        // Delta
        DenseMatrix64F delta = new DenseMatrix64F(dimProcess, dimProcess);
        addTrans(invAS, delta);
        // Variance 2
        CommonOps.addEquals(varianceXX, branchLength, delta);
        // Variance 3
        DenseMatrix64F temp2 = CommonOps.identity(dimProcess);
        CommonOps.addEquals(temp2, -1.0, actualization);
        DenseMatrix64F temp3 = new DenseMatrix64F(dimProcess, dimProcess);
        CommonOps.mult(temp2, inverseSelectionStrength, temp3);
        CommonOps.mult(temp3, delta, temp2);
        addTrans(temp2, temp);
        // All
        CommonOps.addEquals(varianceXX, -1.0, temp);

        //// Variance XY
        DenseMatrix64F varianceXY = new DenseMatrix64F(dimProcess, dimProcess);
        // Variance 1
        CommonOps.multTransB(stationaryVariance, temp3, varianceXY);
        // Variance 2
        CommonOps.mult(temp3, stationaryVariance, temp);
        CommonOps.multTransB(temp, actualization, temp2);
        // All
        CommonOps.addEquals(varianceXY, -1.0, temp2);

        //// Variance YX
        DenseMatrix64F varianceYX = new DenseMatrix64F(dimProcess, dimProcess);
        CommonOps.transpose(varianceXY, varianceYX);

        blockUnwrap(varianceYY, varianceXX, varianceXY, varianceYX, variances, destinationOffset);
        schurComplementInverse(varianceYY, varianceXX, varianceXY, varianceYX, precisions, destinationOffset);
    }

    private void computeIOUActualization(final int scaledOffset,
                                         DenseMatrix64F inverseSelectionStrength) {
        // YY
        DenseMatrix64F actualizationOU = wrap(actualizations, scaledOffset, dimProcess, dimProcess);

        // XX
        DenseMatrix64F temp = CommonOps.identity(dimProcess);
        CommonOps.addEquals(temp, -1.0, actualizationOU);
        DenseMatrix64F actualizationIOU = new DenseMatrix64F(dimProcess, dimProcess);
        CommonOps.mult(inverseSelectionStrength, temp, actualizationIOU);

        // YX and XX
        DenseMatrix64F actualizationYX = new DenseMatrix64F(dimProcess, dimProcess); // zeros
        DenseMatrix64F actualizationXX = CommonOps.identity(dimProcess);

        blockUnwrap(actualizationOU, actualizationXX, actualizationIOU, actualizationYX, actualizations, scaledOffset);
    }

    private void addTrans(DenseMatrix64F A, DenseMatrix64F B) {
        CommonOps.transpose(A, B);
        CommonOps.addEquals(B, A);
    }

    private void blockUnwrap(final DenseMatrix64F YY, final DenseMatrix64F XX, final DenseMatrix64F XY, final DenseMatrix64F YX, final double[] destination, final int offset) {
        for (int i = 0; i < dimProcess; i++) { // Rows
            for (int j = 0; j < dimProcess; j++) {
                destination[offset + i * dimTrait + j] = YY.get(i, j);
                destination[offset + (i + dimProcess) * dimTrait + j + dimProcess] = XX.get(i, j);
            }
            for (int j = 0; j < dimProcess; j++) {
                destination[offset + i * dimTrait + j + dimProcess] = YX.get(i, j);
                destination[offset + (i + dimProcess) * dimTrait + j] = XY.get(i, j);
            }
        }
    }

    private void schurComplementInverse(final DenseMatrix64F A, final DenseMatrix64F D,
                                        final DenseMatrix64F C, final DenseMatrix64F B,
                                        final double[] destination, final int offset) {
        DenseMatrix64F invA = new DenseMatrix64F(dimProcess, dimProcess);
        CommonOps.invert(A, invA);
        DenseMatrix64F invMatD = getSchurInverseComplement(invA, D, C, B);

        DenseMatrix64F invAB = new DenseMatrix64F(dimProcess, dimProcess);
        CommonOps.mult(invA, B, invAB);
        DenseMatrix64F invMatB = new DenseMatrix64F(dimProcess, dimProcess);
        CommonOps.mult(-1.0, invAB, invMatD, invMatB);

        DenseMatrix64F CinvA = new DenseMatrix64F(dimProcess, dimProcess);
        CommonOps.mult(C, invA, CinvA);
        DenseMatrix64F invMatC = new DenseMatrix64F(dimProcess, dimProcess);
        CommonOps.mult(-1.0, invMatD, CinvA, invMatC);

        DenseMatrix64F invMatA = new DenseMatrix64F(dimProcess, dimProcess);
        CommonOps.mult(-1.0, invMatB, CinvA, invMatA);
        CommonOps.addEquals(invMatA, invA);

        blockUnwrap(invMatA, invMatD, invMatC, invMatB, destination, offset);
    }

    private DenseMatrix64F getSchurInverseComplement(final DenseMatrix64F invA, final DenseMatrix64F D,
                                        final DenseMatrix64F C, final DenseMatrix64F B) {
        DenseMatrix64F complement = new DenseMatrix64F(dimProcess, dimProcess);
        DenseMatrix64F tmp = new DenseMatrix64F(dimProcess, dimProcess);
        CommonOps.mult(invA, B, tmp);
        CommonOps.mult(-1.0, C, tmp, complement);
        CommonOps.addEquals(complement, D);
        CommonOps.invert(complement);
        return complement;
    }

    @Override
    public void updateIntegratedOrnsteinUhlenbeckDiffusionMatrices(int precisionIndex, final int[] probabilityIndices,
                                                                   final double[] edgeLengths, final double[] optimalRates,
                                                                   final double[] diagonalStrengthOfSelectionMatrix,
                                                                   final double[] rotation,
                                                                   int updateCount) {

        updateOrnsteinUhlenbeckDiffusionMatrices(precisionIndex, probabilityIndices, edgeLengths, optimalRates,
                diagonalStrengthOfSelectionMatrix, rotation, updateCount);

        if (DEBUG) {
            System.err.println("Matrices (safe with actualized drift, integrated):");
        }

        int matrixTraitSize = dimTrait * dimTrait;
        int matrixProcessSize = dimProcess * dimProcess;
        final int unscaledOffset = matrixProcessSize * precisionIndex;

        final DenseMatrix64F inverseSelectionStrength = getInverseSelectionStrength(diagonalStrengthOfSelectionMatrix, rotation);

        if (TIMING) {
            startTime("drift2");
        }

        int offset = 0;
        for (int up = 0; up < updateCount; ++up) {

            final int pio = dimTrait * probabilityIndices[up];
            final double edgeLength = edgeLengths[up];

            computeIOUActualizedDisplacement(optimalRates, offset, pio, edgeLength, inverseSelectionStrength);
            offset += dimProcess;
        }

        if (TIMING) {
            endTime("drift2");
        }

        if (TIMING) {
            startTime("diffusion2");
        }

        for (int up = 0; up < updateCount; ++up) {

            final int scaledOffset = matrixTraitSize * probabilityIndices[up];
            final double edgeLength = edgeLengths[up];

            computeIOUVarianceBranch(unscaledOffset, scaledOffset, edgeLength, inverseSelectionStrength);
        }

        if (TIMING) {
            endTime("diffusion2");
        }

        if (TIMING) {
            startTime("actualization2");
        }

        for (int up = 0; up < updateCount; ++up) {

            final int scaledOffset = matrixTraitSize * probabilityIndices[up];

            computeIOUActualization(scaledOffset, inverseSelectionStrength);
        }

        if (TIMING) {
            endTime("actualization");
        }

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

/*
 * SafeMultivariateActualizedWithDriftIntegrator.java
 *
 * General OU integrator:
 *   X_child = Q(t) X_parent + displacement(t),
 *   Var[X_child] = Σ - Q Σ Qᵀ
 *
 * Q(t) is a full matrix in trait space. A number of methods are written
 * so that specialisations (e.g. purely diagonal Q) can override only the
 * hooks they need.
 *
 * Copyright © 2002-2024 the BEAST Development Team
 */

package dr.evomodel.treedatalikelihood.continuous.cdi;

import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

import static dr.math.matrixAlgebra.missingData.MissingOps.*;

/**
 * @author Marc A. Suchard
 * @author Paul Bastide
 */
public class SafeMultivariateActualizedWithDriftIntegrator
        extends SafeMultivariateWithDriftIntegrator {

    protected static final boolean DEBUG  = false;
    protected static final boolean TIMING = false;

    /** Stationary variance matrices Σ, one per precision index. */
    protected double[] stationaryVariances;  // dimProcess x dimProcess x diffusionCount

    /** Actualisation matrices Q(t) in trait space, one per branch buffer. */
    protected double[] actualizations;       // dimTrait x dimTrait x bufferCount

    /** Workspaces for tree-traversal operations. */
    protected DenseMatrix64F matrixQdiPip;
    protected DenseMatrix64F matrixQdjPjp;
    protected DenseMatrix64F matrixNiacc;

    /** Whether the eigen-basis rotation is orthonormal. */
    protected final boolean isActualizationSymmetric;

    public SafeMultivariateActualizedWithDriftIntegrator(
            PrecisionType precisionType,
            int numTraits,
            int dimTrait,
            int dimProcess,
            int bufferCount,
            int diffusionCount,
            boolean isActualizationSymmetric) {

        super(precisionType, numTraits, dimTrait, dimProcess, bufferCount, diffusionCount);

        this.isActualizationSymmetric = isActualizationSymmetric;

        allocateStorage();
        System.err.println("Trying SafeMultivariateActualizedWithDriftIntegrator");
    }

    private void allocateStorage() {
        stationaryVariances = new double[dimProcess * dimProcess * diffusionCount];
        actualizations      = new double[dimTrait   * dimTrait   * bufferCount];

        matrixQdiPip = new DenseMatrix64F(dimTrait, dimTrait);
        matrixQdjPjp = new DenseMatrix64F(dimTrait, dimTrait);
        matrixNiacc  = new DenseMatrix64F(dimTrait, 1);
    }

    /* **********************************************************************
     * Branch accessors
     * ******************************************************************* */

    @Override
    public void getBranchActualization(int bufferIndex, double[] actualization) {
        if (bufferIndex == -1) {
            throw new RuntimeException("Not yet implemented");
        }

        assert actualization != null;
        assert actualization.length >= dimTrait * dimTrait;

        System.arraycopy(actualizations,
                bufferIndex * dimTrait * dimTrait,
                actualization, 0,
                dimTrait * dimTrait);
    }

    @Override
    public void getBranch1mActualization(int bufferIndex, double[] actualization) {

        getBranchActualization(bufferIndex, actualization);

        // actualization <- I - actualization
        for (int i = 0; i < dimTrait; i++) {
            for (int j = 0; j < dimTrait; j++) {
                actualization[i * dimTrait + j] = -actualization[i * dimTrait + j];
            }
            actualization[i * dimTrait + i] += 1.0;
        }
    }

    @Override
    public void getBranchExpectation(double[] actualization,
                                     double[] parentValue,
                                     double[] displacement,
                                     double[] expectation) {

        assert expectation   != null && expectation.length   >= dimTrait;
        assert actualization != null && actualization.length >= dimTrait * dimTrait;
        assert parentValue   != null && parentValue.length   >= dimTrait;
        assert displacement  != null && displacement.length  >= dimTrait;

        DenseMatrix64F branchExpectationMatrix = new DenseMatrix64F(dimTrait, 1);
        CommonOps.mult(wrap(actualization, 0, dimTrait, dimTrait),
                wrap(parentValue,   0, dimTrait, 1),
                branchExpectationMatrix);
        CommonOps.addEquals(branchExpectationMatrix,
                wrap(displacement, 0, dimTrait, 1));

        unwrap(branchExpectationMatrix, expectation, 0);
    }

    public double[] getStationaryVariance(int precisionIndex) {
        assert stationaryVariances != null;
        return getMatrixProcess(precisionIndex, stationaryVariances);
    }

//        @Override // TODO add to ContinuousDiffusionIntegrator?
    public void getDiffusionStationaryVariance(int precisionIndex, double[] out) {
        // Σ_stat is stored in stationaryVariances in the ORIGINAL process basis
        double[] src = getStationaryVariance(precisionIndex);
        final int len = dimProcess * dimProcess;
        if (out.length < len) {
            throw new IllegalArgumentException("Output buffer too small for stationary variance");
        }
        System.arraycopy(src, 0, out, 0, len);
    }

//        @Override // TODO add to ContinuousDiffusionIntegrator?
    public void getBranchVariance(int bufferIndex, double[] out) {
        // V is stored in variances as dimTrait × dimTrait per buffer
        final int matrixTraitSize = dimTrait * dimTrait;
        if (out.length < matrixTraitSize) {
            throw new IllegalArgumentException("Output buffer too small for branch variance");
        }
        int offset = matrixTraitSize * bufferIndex;
        System.arraycopy(variances, offset, out, 0, matrixTraitSize);
    }




    /* **********************************************************************
     * Stationary variance Σ
     * ******************************************************************* */

    @Override
    public void setDiffusionStationaryVariance(int precisionIndex,
                                               final double[] alphaEig,
                                               final double[] alphaRot) {

        final int dim = alphaEig.length;
        assert dim == dimProcess;
        assert alphaRot.length == dim * dim;

        final int matrixSize = dim * dim;
        final int offset     = matrixSize * precisionIndex;

        double[] scales = new double[matrixSize];
        scalingMatrix(alphaEig, scales);

        // Fill stationary variance (currently in eigenbasis) via hook.
        setStationaryVariance(offset, scales, matrixSize, alphaRot);

        // Transform back to original space.
        transformMatrixBack(stationaryVariances, offset, alphaRot, 0);

        if (DEBUG) {
            System.err.println("At precision index: " + precisionIndex);
            System.err.println("stationary variance: " +
                    wrap(stationaryVariances, offset, dim, dim));
        }
    }

    /**
     * Hook used by setDiffusionStationaryVariance. The default implementation
     * corresponds to the general (non-diagonal) OU case.
     */

    void setStationaryVariance(int offset,
                               double[] scales,
                               int matrixSize,
                               double[] rotation) {
        assert rotation.length == matrixSize;

        DenseMatrix64F rotMat   = wrap(rotation, 0, dimProcess, dimProcess);
        DenseMatrix64F variance = wrap(inverseDiffusions, offset, dimProcess, dimProcess);

        // Move precision/diffusion into eigenbasis of S.
        transformMatrix(variance, rotMat, isActualizationSymmetric);

        double[] transVar = new double[matrixSize];
        unwrap(variance, transVar, 0);

        // Σ_eig = invDiffusion_eig / (λ_i + λ_j)
        scaleInv(transVar, 0, scales, stationaryVariances, offset, matrixSize);
    }

    /** dest[i] = src[i] / scales[i]. */
    protected static void scaleInv(final double[] source,
                                   final int sourceOffset,
                                   final double[] scales,
                                   final double[] destination,
                                   final int destinationOffset,
                                   final int length) {
        for (int i = 0; i < length; ++i) {
            destination[destinationOffset + i] =
                    source[sourceOffset + i] / scales[i];
        }
    }

    /** scales[i,j] = λ_i + λ_j. */
    protected static void scalingMatrix(double[] eigAlpha, double[] scales) {
        int nEig = eigAlpha.length;
        for (int i = 0; i < nEig; ++i) {
            for (int j = 0; j < nEig; ++j) {
                scales[i * nEig + j] = eigAlpha[i] + eigAlpha[j];
            }
        }
    }

    /* **********************************************************************
     * OU update (non-integrated)
     * ******************************************************************* */

    @Override
    public void updateOrnsteinUhlenbeckDiffusionMatrices(int precisionIndex,
                                                         final int[] probabilityIndices,
                                                         final double[] edgeLengths,
                                                         final double[] optimalRates,
                                                         final double[] diagonalStrengthOfSelectionMatrix,
                                                         final double[] rotation,
                                                         int updateCount) {

        assert diffusions != null;
        assert probabilityIndices.length >= updateCount;
        assert edgeLengths.length        >= updateCount;

        // Let the parent set up diffusions etc.
        super.updateOrnsteinUhlenbeckDiffusionMatrices(precisionIndex, probabilityIndices,
                edgeLengths, optimalRates,
                diagonalStrengthOfSelectionMatrix, rotation, updateCount);

        final int matrixTraitSize   = dimTrait   * dimTrait;
        final int matrixProcessSize = dimProcess * dimProcess;
        final int unscaledOffset    = matrixProcessSize * precisionIndex;

        if (TIMING) startTime("actualization");

        // 1. Q(t) and 1 - Q(t)
        for (int up = 0; up < updateCount; ++up) {
            final double edgeLength     = edgeLengths[up];
            final int   scaledOffsetDiag= dimTrait * probabilityIndices[up];
            final int   scaledOffset    = dimTrait * scaledOffsetDiag;

            computeOUActualization(diagonalStrengthOfSelectionMatrix,
                    rotation,
                    edgeLength,
                    scaledOffsetDiag,
                    scaledOffset);
        }

        if (TIMING) endTime("actualization");
        if (TIMING) startTime("diffusion");

        // 2. Var
        for (int up = 0; up < updateCount; ++up) {
            final double edgeLength     = edgeLengths[up];
            final int   scaledOffset    = matrixTraitSize * probabilityIndices[up];
            final int   scaledOffsetDiag= dimTrait       * probabilityIndices[up];

            computeOUVarianceBranch(unscaledOffset,
                    scaledOffset,
                    scaledOffsetDiag,
                    edgeLength);

            invertVectorSymmPosDef(variances, precisions, scaledOffset, dimProcess);
        }

        if (TIMING) endTime("diffusion");

        assert optimalRates != null;
        assert displacements != null;
        assert optimalRates.length >= updateCount * dimProcess;

        if (TIMING) startTime("drift1");

        // 3. Displacement
        int offset = 0;
        for (int up = 0; up < updateCount; ++up) {
            final int pio          = dimTrait * probabilityIndices[up];
            final int scaledOffset = matrixTraitSize * probabilityIndices[up];

            computeOUActualizedDisplacement(optimalRates,
                    offset,
                    scaledOffset,
                    pio);
            offset += dimProcess;
        }

        if (TIMING) endTime("drift1");
    }

    /**
     * Hook: build Q(t) from diagonal eigenvalues and rotation matrix.
     * Default is full-matrix case.
     */
    void computeOUActualization(final double[] diagonalStrengthOfSelectionMatrix,
                                final double[] rotation,
                                final double edgeLength,
                                final int scaledOffsetDiagonal,
                                final int scaledOffset) {
        double[] diagonalActualizations = new double[dimTrait];

        computeOUDiagonal1mActualization(diagonalStrengthOfSelectionMatrix,
                edgeLength,
                dimTrait,
                diagonalActualizations,
                0);
        oneMinus(diagonalActualizations);

        transformDiagonalMatrixBack(diagonalActualizations,
                actualizations,
                scaledOffset,
                rotation,
                0);
    }

    protected static void computeOUDiagonal1mActualization(final double[] source,
                                                           final double edgeLength,
                                                           final int dim,
                                                           final double[] destination,
                                                           final int destinationOffset) {
        for (int i = 0; i < dim; ++i) {
            destination[destinationOffset + i] =
                    -Math.expm1(-source[i] * edgeLength);
        }
    }

    protected static void oneMinus(double[] x) {
        for (int i = 0; i < x.length; ++i) {
            x[i] = 1.0 - x[i];
        }
    }

    /**
     * Hook: Var[X_child] = Σ - Q Σ Qᵀ (default general case).
     */
    void computeOUVarianceBranch(final int sourceOffset,
                                 final int destinationOffset,
                                 final int destinationOffsetDiagonal,
                                 final double edgeLength) {

        DenseMatrix64F A = wrap(actualizations,      destinationOffset, dimProcess, dimProcess);
        DenseMatrix64F V = wrap(stationaryVariances, sourceOffset,      dimProcess, dimProcess);
        DenseMatrix64F temp = new DenseMatrix64F(dimProcess, dimProcess);

        CommonOps.multTransB(V, A, temp);
        CommonOps.multAdd(-1.0, A, temp, V);

        unwrap(V, variances, destinationOffset);
    }

    /**
     * Hook: (I - Q(t)) μ* (default general case).
     */
    void computeOUActualizedDisplacement(final double[] optimalRates,
                                         final int offset,
                                         final int actualizationOffset,
                                         final int pio) {

        DenseMatrix64F A    = wrap(actualizations, actualizationOffset, dimProcess, dimProcess);
        DenseMatrix64F opt  = wrap(optimalRates,   offset,             dimProcess, 1);
        DenseMatrix64F temp = CommonOps.identity(dimProcess);
        DenseMatrix64F disp = new DenseMatrix64F(dimProcess, 1);

        CommonOps.addEquals(temp, -1.0, A);
        CommonOps.mult(temp, opt, disp);

        unwrap(disp, displacements, pio);
    }

    private static void invertVectorSymmPosDef(final double[] source,
                                               final double[] destination,
                                               final int offset,
                                               final int dim) {
        DenseMatrix64F sourceMatrix      = wrap(source,      offset, dim, dim);
        DenseMatrix64F destinationMatrix = new DenseMatrix64F(dim, dim);

        symmPosDefInvert(sourceMatrix, destinationMatrix);
        unwrap(destinationMatrix, destination, offset);
    }

    /* **********************************************************************
     * Basis transforms (used for stationary variance and integrated OU)
     * ******************************************************************* */

    public static void transformMatrix(DenseMatrix64F matrix,
                                       DenseMatrix64F rotation,
                                       boolean isSymmetric) {
        if (isSymmetric) {
            transformMatrixSymmetric(matrix, rotation);
        } else {
            transformMatrixGeneral(matrix, rotation);
        }
    }

    private static void transformMatrixGeneral(DenseMatrix64F matrix,
                                               DenseMatrix64F rotation) {
        int dim = matrix.getNumRows();
        DenseMatrix64F tmp             = new DenseMatrix64F(dim, dim);
        DenseMatrix64F rotationInverse = new DenseMatrix64F(dim, dim);

        CommonOps.invert(rotation, rotationInverse);
        CommonOps.mult(rotationInverse, matrix, tmp);
        CommonOps.multTransB(tmp, rotationInverse, matrix);
    }

    private static void transformMatrixSymmetric(DenseMatrix64F matrix,
                                                 DenseMatrix64F rotation) {
        int dim = matrix.getNumRows();
        DenseMatrix64F tmp = new DenseMatrix64F(dim, dim);

        CommonOps.multTransA(rotation, matrix, tmp);
        CommonOps.mult(tmp, rotation, matrix);
    }

    private void transformMatrixBack(double[] matrixDouble, int matrixOffset,
                                     double[] rotationDouble, int rotationOffset) {
        DenseMatrix64F matrix   = wrap(matrixDouble,   matrixOffset, dimProcess, dimProcess);
        DenseMatrix64F rotation = wrap(rotationDouble, rotationOffset, dimProcess, dimProcess);
        transformMatrixBack(matrix, rotation);
        unwrap(matrix, matrixDouble, matrixOffset);
    }

    public static void transformMatrixBack(DenseMatrix64F matrix,
                                           DenseMatrix64F rotation) {
        int dim = matrix.getNumRows();
        DenseMatrix64F tmp = new DenseMatrix64F(dim, dim);

        CommonOps.multTransB(matrix, rotation, tmp);
        CommonOps.mult(rotation, tmp, matrix);
    }

    private void transformDiagonalMatrixBack(double[] diagonalMatrix,
                                             double[] matrixDestination,
                                             int matrixOffset,
                                             double[] rotationDouble,
                                             int rotationOffset) {
        DenseMatrix64F matrix   = wrapDiagonal(diagonalMatrix, matrixOffset, dimProcess);
        DenseMatrix64F rotation = wrap(rotationDouble, rotationOffset, dimProcess, dimProcess);
        transformMatrixBase(matrix, rotation);
        unwrap(matrix, matrixDestination, matrixOffset);
    }

    private DenseMatrix64F getInverseSelectionStrength(double[] diagonalMatrix,
                                                       double[] rotationDouble) {
        DenseMatrix64F matrix   = wrapDiagonalInverse(diagonalMatrix, 0, dimProcess);
        DenseMatrix64F rotation = wrap(rotationDouble, 0, dimProcess, dimProcess);
        transformMatrixBase(matrix, rotation);
        return matrix;
    }

    void transformMatrixBase(DenseMatrix64F matrix,
                             DenseMatrix64F rotation) {
        if (isActualizationSymmetric) {
            transformMatrixBack(matrix, rotation);
        } else {
            transformMatrixBaseGeneral(matrix, rotation);
        }
    }

    private void transformMatrixBaseGeneral(DenseMatrix64F matrix,
                                            DenseMatrix64F rotation) {
        DenseMatrix64F tmp = new DenseMatrix64F(dimProcess, dimProcess);
        CommonOps.mult(rotation, matrix, tmp);
        CommonOps.invert(rotation); // Warning: side effect on rotation.
        CommonOps.mult(tmp, rotation, matrix);
    }

    /* **********************************************************************
     * Integrated OU (same as in your original general class)
     * ******************************************************************* */

    @Override
    public void updateIntegratedOrnsteinUhlenbeckDiffusionMatrices(
            int precisionIndex,
            final int[] probabilityIndices,
            final double[] edgeLengths,
            final double[] optimalRates,
            final double[] diagonalStrengthOfSelectionMatrix,
            final double[] rotation,
            int updateCount) {

        super.updateOrnsteinUhlenbeckDiffusionMatrices(precisionIndex, probabilityIndices,
                edgeLengths, optimalRates,
                diagonalStrengthOfSelectionMatrix, rotation, updateCount);

        int matrixTraitSize   = dimTrait   * dimTrait;
        int matrixProcessSize = dimProcess * dimProcess;
        final int unscaledOffset = matrixProcessSize * precisionIndex;

        final DenseMatrix64F inverseSelectionStrength =
                getInverseSelectionStrength(diagonalStrengthOfSelectionMatrix, rotation);

        if (TIMING) startTime("drift2");

        int offset = 0;
        for (int up = 0; up < updateCount; ++up) {

            final int    pio        = dimTrait * probabilityIndices[up];
            final double edgeLength = edgeLengths[up];

            computeIOUActualizedDisplacement(optimalRates, offset, pio,
                    edgeLength, inverseSelectionStrength);
            offset += dimProcess;
        }

        if (TIMING) endTime("drift2");
        if (TIMING) startTime("diffusion2");

        for (int up = 0; up < updateCount; ++up) {

            final int    scaledOffset = matrixTraitSize * probabilityIndices[up];
            final double edgeLength   = edgeLengths[up];

            computeIOUVarianceBranch(unscaledOffset, scaledOffset,
                    edgeLength, inverseSelectionStrength);
        }

        if (TIMING) endTime("diffusion2");
        if (TIMING) startTime("actualization2");

        for (int up = 0; up < updateCount; ++up) {
            final int scaledOffset = matrixTraitSize * probabilityIndices[up];
            computeIOUActualization(scaledOffset, inverseSelectionStrength);
        }

        if (TIMING) endTime("actualization2");
    }

    void computeIOUActualizedDisplacement(final double[] optimalRates,
                                          final int offset,
                                          final int pio,
                                          double branchLength,
                                          DenseMatrix64F inverseSelectionStrength) {
        DenseMatrix64F displacementOU = wrap(displacements, pio, dimProcess, 1);
        DenseMatrix64F optVal         = wrap(optimalRates, offset, dimProcess, 1);
        DenseMatrix64F displacement   = new DenseMatrix64F(dimProcess, 1);

        CommonOps.mult(inverseSelectionStrength, displacementOU, displacement);
        CommonOps.scale(-1.0, displacement);
        CommonOps.addEquals(displacement, branchLength, optVal);

        unwrap(displacement, displacements, pio + dimProcess);
    }

    void computeIOUVarianceBranch(final int sourceOffset,
                                  final int destinationOffset,
                                  double branchLength,
                                  DenseMatrix64F inverseSelectionStrength) {

        DenseMatrix64F actualization      = wrap(actualizations,      destinationOffset, dimProcess, dimProcess);
        DenseMatrix64F stationaryVariance = wrap(stationaryVariances, sourceOffset,      dimProcess, dimProcess);

        DenseMatrix64F invAS = new DenseMatrix64F(dimProcess, dimProcess);
        CommonOps.mult(inverseSelectionStrength, stationaryVariance, invAS);

        // Var(YY)
        DenseMatrix64F varianceYY = wrap(variances, destinationOffset, dimProcess, dimProcess);

        // Var(XX)
        DenseMatrix64F varianceXX = new DenseMatrix64F(dimProcess, dimProcess);
        CommonOps.multTransB(invAS, inverseSelectionStrength, varianceXX);

        DenseMatrix64F temp = new DenseMatrix64F(dimProcess, dimProcess);
        CommonOps.multTransB(varianceXX, actualization, temp);
        CommonOps.multAdd(-1.0, actualization, temp, varianceXX);

        DenseMatrix64F delta = new DenseMatrix64F(dimProcess, dimProcess);
        addTrans(invAS, delta);
        CommonOps.addEquals(varianceXX, branchLength, delta);

        DenseMatrix64F temp2 = CommonOps.identity(dimProcess);
        CommonOps.addEquals(temp2, -1.0, actualization);
        DenseMatrix64F temp3 = new DenseMatrix64F(dimProcess, dimProcess);
        CommonOps.mult(temp2, inverseSelectionStrength, temp3);
        CommonOps.mult(temp3, delta, temp2);
        addTrans(temp2, temp);
        CommonOps.addEquals(varianceXX, -1.0, temp);

        // Var(XY)
        DenseMatrix64F varianceXY = new DenseMatrix64F(dimProcess, dimProcess);
        CommonOps.multTransB(stationaryVariance, temp3, varianceXY);

        CommonOps.mult(temp3, stationaryVariance, temp);
        CommonOps.multTransB(temp, actualization, temp2);
        CommonOps.addEquals(varianceXY, -1.0, temp2);

        // Var(YX)
        DenseMatrix64F varianceYX = new DenseMatrix64F(dimProcess, dimProcess);
        CommonOps.transpose(varianceXY, varianceYX);

        blockUnwrap(varianceYY, varianceXX, varianceXY, varianceYX,
                variances, destinationOffset);
        schurComplementInverse(varianceYY, varianceXX, varianceXY, varianceYX,
                precisions, destinationOffset);
    }

    void computeIOUActualization(final int scaledOffset,
                                 DenseMatrix64F inverseSelectionStrength) {
        DenseMatrix64F actualizationOU = wrap(actualizations, scaledOffset, dimProcess, dimProcess);

        DenseMatrix64F temp = CommonOps.identity(dimProcess);
        CommonOps.addEquals(temp, -1.0, actualizationOU);
        DenseMatrix64F actualizationIOU = new DenseMatrix64F(dimProcess, dimProcess);
        CommonOps.mult(inverseSelectionStrength, temp, actualizationIOU);

        DenseMatrix64F actualizationYX = new DenseMatrix64F(dimProcess, dimProcess); // zeros
        DenseMatrix64F actualizationXX = CommonOps.identity(dimProcess);

        blockUnwrap(actualizationOU, actualizationXX, actualizationIOU, actualizationYX,
                actualizations, scaledOffset);
    }

    private void addTrans(DenseMatrix64F A, DenseMatrix64F B) {
        CommonOps.transpose(A, B);
        CommonOps.addEquals(B, A);
    }

    private void blockUnwrap(final DenseMatrix64F YY,
                             final DenseMatrix64F XX,
                             final DenseMatrix64F XY,
                             final DenseMatrix64F YX,
                             final double[] destination,
                             final int offset) {

        for (int i = 0; i < dimProcess; i++) {
            for (int j = 0; j < dimProcess; j++) {
                destination[offset + i * dimTrait + j] =
                        YY.get(i, j);
                destination[offset + (i + dimProcess) * dimTrait + j + dimProcess] =
                        XX.get(i, j);
            }
            for (int j = 0; j < dimProcess; j++) {
                destination[offset + i * dimTrait + j + dimProcess] =
                        YX.get(i, j);
                destination[offset + (i + dimProcess) * dimTrait + j] =
                        XY.get(i, j);
            }
        }
    }

    private void schurComplementInverse(final DenseMatrix64F A,
                                        final DenseMatrix64F D,
                                        final DenseMatrix64F C,
                                        final DenseMatrix64F B,
                                        final double[] destination,
                                        final int offset) {

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

        blockUnwrap(invMatA, invMatD, invMatC, invMatB,
                destination, offset);
    }

    private DenseMatrix64F getSchurInverseComplement(final DenseMatrix64F invA,
                                                     final DenseMatrix64F D,
                                                     final DenseMatrix64F C,
                                                     final DenseMatrix64F B) {
        DenseMatrix64F complement = new DenseMatrix64F(dimProcess, dimProcess);
        DenseMatrix64F tmp        = new DenseMatrix64F(dimProcess, dimProcess);

        CommonOps.mult(invA, B, tmp);
        CommonOps.mult(-1.0, C, tmp, complement);
        CommonOps.addEquals(complement, D);
        CommonOps.invert(complement);
        return complement;
    }

    /* **********************************************************************
     * Tree traversal – general case
     * ******************************************************************* */

    @Override
    void actualizePrecision(DenseMatrix64F Pjp, DenseMatrix64F QjPjp,
                            int jbo, int jmo, int jdo) {
        final DenseMatrix64F Qdj = wrap(actualizations, jmo, dimTrait, dimTrait);
        scalePrecision(Qdj, Pjp, QjPjp, Pjp);
    }

    @Override
    void actualizeVariance(DenseMatrix64F Vip, int ibo, int imo, int ido) {
        final DenseMatrix64F Qdi  = wrap(actualizations, imo, dimTrait, dimTrait);
        final DenseMatrix64F QiVp = matrixQdiPip;
        scaleVariance(Qdi, Vip, QiVp, Vip);
    }

    @Override
    void scaleAndDriftMean(int ibo, int imo, int ido) {
        final DenseMatrix64F Qdi = wrap(actualizations, imo, dimTrait, dimTrait);
        final DenseMatrix64F ni  = wrap(preOrderPartials, ibo, dimTrait, 1);

        CommonOps.mult(Qdi, ni, matrixNiacc);
        unwrap(matrixNiacc, preOrderPartials, ibo);

        for (int g = 0; g < dimTrait; ++g) {
            preOrderPartials[ibo + g] += displacements[ido + g];
        }
    }

    @Override
    void computePartialPrecision(int ido, int jdo, int imo, int jmo,
                                 DenseMatrix64F Pip,
                                 DenseMatrix64F Pjp,
                                 DenseMatrix64F Pk) {

        final DenseMatrix64F Qdi = wrap(actualizations, imo, dimTrait, dimTrait);
        final DenseMatrix64F Qdj = wrap(actualizations, jmo, dimTrait, dimTrait);

        final DenseMatrix64F QdiPip  = matrixQdiPip;
        final DenseMatrix64F QdiPipQ = matrix0;
        scalePrecision(Qdi, Pip, QdiPip, QdiPipQ);

        final DenseMatrix64F QdjPjpQ = matrix1;
        final DenseMatrix64F QdjPjp  = matrixQdjPjp;
        scalePrecision(Qdj, Pjp, QdjPjp, QdjPjpQ);

        CommonOps.add(QdiPipQ, QdjPjpQ, Pk);
    }

    private void scalePrecision(DenseMatrix64F Q,
                                DenseMatrix64F P,
                                DenseMatrix64F QtP,
                                DenseMatrix64F QtPQ) {
        CommonOps.multTransA(Q, P, QtP);
        CommonOps.mult(QtP, Q, QtPQ);
        forceSymmetric(QtPQ);
    }

    private void scaleVariance(DenseMatrix64F Q,
                               DenseMatrix64F V,
                               DenseMatrix64F QV,
                               DenseMatrix64F QVQt) {
        CommonOps.mult(Q, V, QV);
        CommonOps.multTransB(QV, Q, QVQt);
    }

    @Override
    void computeWeightedSum(final double[] ipartial,
                            final double[] jpartial,
                            final int dimTrait,
                            final double[] out) {
        weightedSum(ipartial, 0, matrixQdiPip,
                jpartial, 0, matrixQdjPjp,
                dimTrait, out);
    }
}

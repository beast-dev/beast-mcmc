package dr.evomodel.treedatalikelihood.continuous.cdi;

import dr.math.matrixAlgebra.WrappedVector;
import dr.math.matrixAlgebra.missingData.InversionResult;

import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

import static dr.math.matrixAlgebra.missingData.InversionResult.Code.NOT_OBSERVED;
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
    public void getBranchMatrices(int bufferIndex, double[] precision, double[] displacement, double[] actualization) {
        if (bufferIndex == -1) {
            throw new RuntimeException("Not yet implemented");
        }

        assert (precision != null);
        assert (precision.length >= dimTrait * dimTrait);

        assert (displacement != null);
        assert (displacement.length >= dimTrait);

        assert (actualization != null);
        assert (actualization.length >= dimTrait * dimTrait);

        System.arraycopy(precisions, bufferIndex * dimTrait * dimTrait,
                precision, 0, dimTrait * dimTrait);

        System.arraycopy(displacements, bufferIndex * dimTrait,
                displacement, 0, dimTrait);

        System.arraycopy(actualizations, bufferIndex * dimTrait * dimTrait,
                actualization, 0, dimTrait * dimTrait);
    }

    private static final boolean TIMING = false;

    private double[] vector1;
    private double[] vector2;

    private void allocateStorage() {

        displacements = new double[dimTrait * bufferCount];
        precisions = new double[dimTrait * dimTrait * bufferCount];
        variances = new double[dimTrait * dimTrait * bufferCount];
        vector1 = new double[dimTrait];
        vector2 = new double[dimTrait];

        actualizations = new double[dimTrait * dimTrait * bufferCount];
        stationaryVariances = new double[dimTrait * dimTrait * diffusionCount];

        matrix7 = new DenseMatrix64F(dimTrait, dimTrait);
        matrix8 = new DenseMatrix64F(dimTrait, dimTrait);
        matrix9 = new DenseMatrix64F(dimTrait, dimTrait);
        matrix10 = new DenseMatrix64F(dimTrait, dimTrait);
    }

    @Override
    public void setDiffusionStationaryVariance(int precisionIndex, final double[] alphaEig, final double[] alphaRot) {

        super.setDiffusionStationaryVariance(precisionIndex, alphaEig, alphaRot);

        // Transform back in original space
        final int offset = dimTrait * dimTrait * precisionIndex;
        DenseMatrix64F stationaryVariance = wrap(stationaryVariances, offset, dimTrait, dimTrait);
        DenseMatrix64F rot = wrap(alphaRot, 0, dimTrait, dimTrait);
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

    @Override
    protected void updatePartial(
            final int kBuffer,
            final int iBuffer,
            final int iMatrix,
            final int jBuffer,
            final int jMatrix,
            final boolean incrementOuterProducts
    ) {

        if (incrementOuterProducts) {
            throw new RuntimeException("Outer-products are not supported.");
        }

        if (TIMING) {
            startTime("total");
        }

        // Determine buffer offsets
        int kbo = dimPartial * kBuffer;
        int ibo = dimPartial * iBuffer;
        int jbo = dimPartial * jBuffer;

        // Determine matrix offsets
        final int imo = dimTrait * dimTrait * iMatrix;
        final int jmo = dimTrait * dimTrait * jMatrix;

        // Read variance increments along descendent branches of k

        // TODO Fix
//        final double vi = variances[imo];
//        final double vj = variances[jmo];

        final DenseMatrix64F Vd = wrap(inverseDiffusions, precisionOffset, dimTrait, dimTrait);
//        final DenseMatrix64F Pd = wrap(diffusions, precisionOffset, dimTrait, dimTrait);

        final DenseMatrix64F Vdi = wrap(variances, imo, dimTrait, dimTrait);
        final DenseMatrix64F Vdj = wrap(variances, jmo, dimTrait, dimTrait);

        final DenseMatrix64F Pdi = wrap(precisions, imo, dimTrait, dimTrait); // TODO Only if needed
        final DenseMatrix64F Pdj = wrap(precisions, jmo, dimTrait, dimTrait); // TODO Only if needed

        final DenseMatrix64F Qdi = wrap(actualizations, imo, dimTrait, dimTrait);
        final DenseMatrix64F Qdj = wrap(actualizations, jmo, dimTrait, dimTrait);

        // TODO End fix

        if (DEBUG) {
            System.err.println("variance diffusion: " + Vd);
//            System.err.println("\tvi: " + vi + " vj: " + vj);
            System.err.println("precisionOffset = " + precisionOffset);
            System.err.println("\tVdi: " + Vdi);
            System.err.println("\tVdj: " + Vdj);
            System.err.println("\tPdi: " + Pdi);
            System.err.println("\tPdj: " + Pdj);
        }

        // For each trait // TODO in parallel
        for (int trait = 0; trait < numTraits; ++trait) {

            // Layout, offset, dim
            // trait, 0, dT
            // precision, dT, dT * dT
            // variance, dT + dT * dT, dT * dT
            // scalar, dT + 2 * dT * dT, 1

            if (TIMING) {
                startTime("peel1");
            }

            // Increase variance along the branches i -> k and j -> k

            // A. Get current precision of i and j
//            final double lpi = partials[ibo + dimTrait + 2 * dimTrait * dimTrait];
//            final double lpj = partials[jbo + dimTrait + 2 * dimTrait * dimTrait];

            final DenseMatrix64F Pi = wrap(partials, ibo + dimTrait, dimTrait, dimTrait);
            final DenseMatrix64F Pj = wrap(partials, jbo + dimTrait, dimTrait, dimTrait);

            if (TIMING) {
                endTime("peel1");
                startTime("peel2");
            }

            // B. Integrate along branch using two matrix inversions
//            final double lpip = Double.isInfinite(lpi) ?
//                    1.0 / vi : lpi / (1.0 + lpi * vi);
//            final double lpjp = Double.isInfinite(lpj) ?
//                    1.0 / vj : lpj / (1.0 + lpj * vj);

            InversionResult ci;
            InversionResult cj;

            final DenseMatrix64F Pip = matrix2;
            final DenseMatrix64F Pjp = matrix3;

//            boolean useVariance = anyDiagonalInfinities(Pi) || anyDiagonalInfinities(Pj);
            final boolean useVariancei = anyDiagonalInfinities(Pi);
            final boolean useVariancej = anyDiagonalInfinities(Pj);

            if (useVariancei) {

                final DenseMatrix64F Vip = matrix0;
                final DenseMatrix64F Vi = wrap(partials, ibo + dimTrait + dimTrait * dimTrait, dimTrait, dimTrait);
//                CommonOps.add(Vi, vi, Vd, Vip);  // TODO Fix
                CommonOps.add(Vi, Vdi, Vip);
                ci = safeInvert(Vip, Pip, true);

            } else {

                final DenseMatrix64F PiPlusPd = matrix0;
//                CommonOps.add(Pi, 1.0 / vi, Pd, PiPlusPd); // TODO Fix
                CommonOps.add(Pi, Pdi, PiPlusPd);
                final DenseMatrix64F PiPlusPdInv = new DenseMatrix64F(dimTrait, dimTrait);
                safeInvert(PiPlusPd, PiPlusPdInv, false);
                CommonOps.mult(PiPlusPdInv, Pi, Pip);
                CommonOps.mult(Pi, Pip, PiPlusPdInv);
                CommonOps.add(Pi, -1, PiPlusPdInv, Pip);
                ci = safeDeterminant(Pip, false);
            }
            // Actualization
            final DenseMatrix64F QdiPip = matrix7;
            CommonOps.multTransA(Qdi, Pip, QdiPip);

            if (useVariancej) {

                final DenseMatrix64F Vjp = matrix1;
                final DenseMatrix64F Vj = wrap(partials, jbo + dimTrait + dimTrait * dimTrait, dimTrait, dimTrait);
//                CommonOps.add(Vj, vj, Vd, Vjp); // TODO Fix
                CommonOps.add(Vj, Vdj, Vjp);
                cj = safeInvert(Vjp, Pjp, true);

            } else {

                final DenseMatrix64F PjPlusPd = matrix1;
//                CommonOps.add(Pj, 1.0 / vj, Pd, PjPlusPd); // TODO Fix
                CommonOps.add(Pj, Pdj, PjPlusPd);
                final DenseMatrix64F PjPlusPdInv = new DenseMatrix64F(dimTrait, dimTrait);
                safeInvert(PjPlusPd, PjPlusPdInv, false);
                CommonOps.mult(PjPlusPdInv, Pj, Pjp);
                CommonOps.mult(Pj, Pjp, PjPlusPdInv);
                CommonOps.add(Pj, -1, PjPlusPdInv, Pjp);
                cj = safeDeterminant(Pjp, false);
            }
            // Actualization
            final DenseMatrix64F QdjPjp = matrix8;
            CommonOps.multTransA(Qdj, Pjp, QdjPjp);

            if (TIMING) {
                endTime("peel2");
                startTime("peel3");
            }

            // Compute partial mean and precision at node k

            // A. Partial precision and variance (for later use) using one matrix inversion
//            final double lpk = lpip + lpjp;

//                final DenseMatrix64F Pk = new DenseMatrix64F(dimTrait, dimTrait);
            final DenseMatrix64F Pk = matrix4;
            final DenseMatrix64F QdiPipQdi = matrix9; // TODO Could re-use matrix0 and matrix1 ?
            final DenseMatrix64F QdjPjpQdj = matrix10;
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

//                final DenseMatrix64F Vk = new DenseMatrix64F(dimTrait, dimTrait);
//            final DenseMatrix64F Vk = matrix5;

//            if (useVariance) {
////            InversionResult ck =
//                    safeInvert(Pk, Vk, true);
//            }
//            InversionResult ck = safeDeterminant(Pk, true);

//            System.err.println(ck);
//            System.err.println(ck2);
//            System.exit(-1);

            // B. Partial mean

            if (TIMING) {
                endTime("peel3");
                startTime("peel4");
            }

            final double[] displacementi = vector1;
            final double[] displacementj = vector2;

            final int ido = dimTrait * iMatrix;
            final int jdo = dimTrait * jMatrix;

            for (int g = 0; g < dimTrait; ++g) {
                displacementi[g] = partials[ibo + g] - displacements[ido + g];
                displacementj[g] = partials[jbo + g] - displacements[jdo + g];
            }

//                final double[] tmp = new double[dimTrait];
            final double[] tmp = vector0;
            for (int g = 0; g < dimTrait; ++g) {
                double sum = 0.0;
                for (int h = 0; h < dimTrait; ++h) {
                    sum += QdiPip.unsafe_get(g, h) * displacementi[h];
                    sum += QdjPjp.unsafe_get(g, h) * displacementj[h];
                }
                tmp[g] = sum;
            }

//            for (int g = 0; g < dimTrait; ++g) {
//                double sum = 0.0;
//                for (int h = 0; h < dimTrait; ++h) {
//                    sum += Vk.unsafe_get(g, h) * tmp[h];
//                }
//                partials[kbo + g] = sum;
//            }

            final WrappedVector kPartials = new WrappedVector.Raw(partials, kbo, dimTrait);
            final WrappedVector wrapTmp = new WrappedVector.Raw(tmp, 0, dimTrait);

//            System.err.println(kPartials);
//            System.err.println(ck.getDeterminant());

            InversionResult ck = safeSolve(Pk, wrapTmp, kPartials, true);

//            System.err.println(kPartials);
//            System.err.println(ck.getDeterminant());
//            System.exit(-1);

            if (TIMING) {
                endTime("peel4");
                startTime("peel5");
            }

            // C. Store precision
//            partials[kbo + dimTrait + 2 * dimTrait * dimTrait] = lpk;

            unwrap(Pk, partials, kbo + dimTrait);

//            if (useVariance) {
//                unwrap(Vk, partials, kbo + dimTrait + dimTrait * dimTrait);
//            }

            if (TIMING) {
                endTime("peel5");
            }

            if (DEBUG) {
                System.err.println("\ttrait: " + trait);
                System.err.println("Pi: " + Pi);
                System.err.println("Pj: " + Pj);
                System.err.println("Pk: " + Pk);
                System.err.print("\t\tmean i:");
                for (int e = 0; e < dimTrait; ++e) {
                    System.err.print(" " + partials[ibo + e]);
                }
                System.err.print("\t\tdisp i:");
                for (int e = 0; e < dimTrait; ++e) {
                    System.err.print(" " + displacements[ido + e]);
                }
                System.err.println("");
                System.err.print("\t\tmean j:");
                for (int e = 0; e < dimTrait; ++e) {
                    System.err.print(" " + partials[jbo + e]);
                }
                System.err.print("\t\tdisp j:");
                for (int e = 0; e < dimTrait; ++e) {
                    System.err.print(" " + displacements[jdo + e]);
                }
                System.err.println("");
                System.err.print("\t\tmean k:");
                for (int e = 0; e < dimTrait; ++e) {
                    System.err.print(" " + partials[kbo + e]);
                }
                System.err.println("");
            }

            // Computer remainder at node k
            double remainder = 0.0;

            if (DEBUG) {
                System.err.println("i status: " + ci);
                System.err.println("j status: " + cj);
                System.err.println("k status: " + ck);
                System.err.println("Pip: " + Pip);
//                System.err.println("Vip: " + Vip);
                System.err.println("Pjp: " + Pjp);
//                System.err.println("Vjp: " + Vjp);
            }

            if (!(ci.getReturnCode() == NOT_OBSERVED || cj.getReturnCode() == NOT_OBSERVED)) {
//                if (ci == InversionReturnCode.FULLY_OBSERVED && cj == InversionReturnCode.FULLY_OBSERVED) {
                // TODO Fix for partially observed
//                if (pi != 0 && pj != 0) {
//
                if (TIMING) {
                    startTime("remain");
                }

                // Inner products
                double SSk = 0;
                double SSj = 0;
                double SSi = 0;

                // vector-matrix-vector TODO in parallel
                for (int g = 0; g < dimTrait; ++g) {
//                    final double ig = partials[ibo + g];
//                    final double jg = partials[jbo + g];
                    final double ig = displacementi[g];
                    final double jg = displacementj[g];
                    final double kg = partials[kbo + g];

                    for (int h = 0; h < dimTrait; ++h) {
//                        final double ih = partials[ibo + h];
//                        final double jh = partials[jbo + h];
                        final double ih = displacementi[h];
                        final double jh = displacementj[h];
                        final double kh = partials[kbo + h];

                        SSi += ig * Pip.unsafe_get(g, h) * ih;
                        SSj += jg * Pjp.unsafe_get(g, h) * jh;
                        SSk += kg * Pk.unsafe_get(g, h) * kh;
                    }
                }
                // TODO Can get SSi + SSj - SSk from inner product w.r.t Pt (see outer-products below)?

                if (DEBUG) {
                    System.err.println("\t\t\tSSi = " + (SSi));
                    System.err.println("\t\t\tSSj = " + (SSj));
                    System.err.println("\t\t\tSSk = " + (SSk));
                }

                remainder += -0.5 * (SSi + SSj - SSk);

            }

            int dimensionChange = ci.getEffectiveDimension() + cj.getEffectiveDimension()
                    - ck.getEffectiveDimension();

            remainder += -dimensionChange * LOG_SQRT_2_PI;

            double deti = 0;
            double detj = 0;
            double detk = 0;
            if (!(ci.getReturnCode() == NOT_OBSERVED)) {
                deti = Math.log(ci.getDeterminant()); // TODO: use det(exp(M)) = exp(tr(M)) ? (Qdi = exp(-A l_i))
            }
            if (!(cj.getReturnCode() == NOT_OBSERVED)) {
                detj = Math.log(cj.getDeterminant());
            }
            if (!(ck.getReturnCode() == NOT_OBSERVED)) {
                detk = Math.log(ck.getDeterminant());
            }
            remainder += - 0.5 * (deti + detj + detk);

            if (DEBUG) {
                System.err.println("\t\t\tdeti = " + Math.log(ci.getDeterminant()));
                System.err.println("\t\t\tdetj = " + Math.log(cj.getDeterminant()));
                System.err.println("\t\t\tdetk = " + Math.log(ck.getDeterminant()));
                System.err.println("\t\tremainder: " + remainder);
//                        System.exit(-1);
            }

            if (TIMING) {
                endTime("remain");
            }

            // Accumulate remainder up tree and store

            remainders[kBuffer * numTraits + trait] = remainder
                    + remainders[iBuffer * numTraits + trait] + remainders[jBuffer * numTraits + trait];

            // Get ready for next trait
            kbo += dimPartialForTrait;
            ibo += dimPartialForTrait;
            jbo += dimPartialForTrait;

        }

        if (TIMING) {
            endTime("total");
        }
    }

    private double[] actualizations;
//    private double[] stationaryVariances;
    private DenseMatrix64F matrix7;
    private DenseMatrix64F matrix8;
    private DenseMatrix64F matrix9;
    private DenseMatrix64F matrix10;

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

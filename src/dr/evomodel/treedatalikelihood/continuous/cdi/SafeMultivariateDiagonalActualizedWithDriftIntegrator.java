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

public class SafeMultivariateDiagonalActualizedWithDriftIntegrator extends SafeMultivariateWithDriftIntegrator {

    private static boolean DEBUG = false;

    public SafeMultivariateDiagonalActualizedWithDriftIntegrator(PrecisionType precisionType, int numTraits, int dimTrait, int bufferCount,
                                                         int diffusionCount) {
        super(precisionType, numTraits, dimTrait, bufferCount, diffusionCount);

        allocateStorage();

        System.err.println("Trying SafeMultivariateDiagonalActualizedWithDriftIntegrator");
    }

    // NOTE TO PB: need to merge all SafeMultivariate* together and then delegate specialized work ... avoid massive code duplication

    @Override
    public void getBranchMatrices(int bufferIndex, double[] precision, double[] displacement, double[] diagonalActualization) {
        if (bufferIndex == -1) {
            throw new RuntimeException("Not yet implemented");
        }

        assert (precision != null);
        assert (precision.length >= dimTrait * dimTrait);

        assert (displacement != null);
        assert (displacement.length >= dimTrait);

        assert (diagonalActualization != null);
        assert (diagonalActualization.length >= dimTrait);

        System.arraycopy(precisions, bufferIndex * dimTrait * dimTrait,
                precision, 0, dimTrait * dimTrait);

        System.arraycopy(displacements, bufferIndex * dimTrait,
                displacement, 0, dimTrait);

        System.arraycopy(diagonalActualizations, bufferIndex * dimTrait,
                diagonalActualization, 0, dimTrait);
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

        diagonalActualizations = new double[dimTrait * bufferCount];
        stationaryVariances = new double[dimTrait * dimTrait * diffusionCount];

//        matrix7 = new DenseMatrix64F(dimTrait, dimTrait);
//        matrix8 = new DenseMatrix64F(dimTrait, dimTrait);
        matrix9 = new DenseMatrix64F(dimTrait, dimTrait);
        matrix10 = new DenseMatrix64F(dimTrait, dimTrait);
    }


    @Override
    public void setDiffusionStationaryVariance(int precisionIndex, final double[] diagonalAlpha, final double[] rotation) {

        assert (stationaryVariances != null);

        final int matrixSize = dimTrait * dimTrait;
        final int offset = matrixSize * precisionIndex;

        double[] scales = new double[matrixSize];
        scalingMatrix(diagonalAlpha, scales);

        if (rotation.length == 0){
            scaleInv(inverseDiffusions, offset, scales, stationaryVariances, offset, matrixSize);
        } else {
            assert(rotation.length == matrixSize);
            DenseMatrix64F rotMat = wrap(rotation, 0, dimTrait, dimTrait);
            DenseMatrix64F variance = wrap(inverseDiffusions, offset, dimTrait, dimTrait);
            transformMatrix(variance, rotMat);
            double[] transVar = new double[matrixSize];
            unwrap(variance, transVar, 0);
            scaleInv(transVar, 0, scales, stationaryVariances, offset, matrixSize);
        }
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

    private static void scale(final double[] source,
                              final int sourceOffset,
                              final double[] scales,
                              final double[] destination,
                              final int destinationOffset,
                              final int length) {
        for (int i = 0; i < length; ++i) {
            destination[destinationOffset + i] = scales[i] * source[sourceOffset + i];
        }
    }

    private static void scaleInv(final double[] source,
                                 final int sourceOffset,
                                 final double[] scales,
                                 final double[] destination,
                                 final int destinationOffset,
                                 final int length) {
        for (int i = 0; i < length; ++i) {
            destination[destinationOffset + i] =  source[sourceOffset + i] / scales[i];
        }
    }

    private static void scalingMatrix(double[] eigAlpha, double[] scales){
        int nEig = eigAlpha.length;
        for (int i = 0; i < nEig; ++i) {
            for (int j = 0; j < nEig; ++j) {
                scales[i * nEig + j] = eigAlpha[i] + eigAlpha[j];
            }
        }
    }

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

            final int scaledOffset = dimTrait * probabilityIndices[up];

            computeDiagonalActualization(diagonalStrengthOfSelectionMatrix, edgeLength, dimTrait, diagonalActualizations, scaledOffset);
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
            final int scaledOffsetDiagonal = dimTrait * probabilityIndices[up];

            double[] scales = new double[matrixSize];
            scalingActualizationMatrix(diagonalActualizations, scaledOffsetDiagonal, dimTrait, scales);

            scale(stationaryVariances, unscaledOffset, scales, variances, scaledOffset, matrixSize);
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

//            final int scaledOffset = matrixSize * probabilityIndices[up];
            final int pio = dimTrait * probabilityIndices[up];

            double[] scales = new double[dimTrait];
            scalingActualizationDisplacement(diagonalActualizations, pio, dimTrait, scales);

            scale(optimalRates, offset, scales, displacements, pio, dimTrait);
            offset += dimTrait;
        }

        if (TIMING) {
            endTime("drift1");
        }

        // NOTE TO PB: very complex function, why multiple for (up = 0; up < updateCount; ++up) ?

        precisionOffset = dimTrait * dimTrait * precisionIndex;
        precisionLogDet = determinants[precisionIndex];
    }

    protected static void computeDiagonalActualization(final double[] source,
                                                     final double edgeLength,
                                                     final int dim,
                                                     final double[] destination,
                                                     final int destinationOffset) {
        for (int i = 0; i < dim; ++i) {
            destination[destinationOffset + i] = Math.exp(- source[i] * edgeLength);
        }
    }

    private static void scalingActualizationMatrix(final double[] source,
                                                   final int offset,
                                                   final int dim,
                                                   final double[] destination) {
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                destination[i * dim + j] = 1 - source[offset + i] * source[offset + j];
            }
        }
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


    private static void scalingActualizationDisplacement(final double[] source,
                                                         final int offset,
                                                         final int dim,
                                                         final double[] destination) {
        for (int i = 0; i < dim; ++i) {
                destination[i] = 1 - source[offset + i];
        }
    }

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

        // Determine diagonal matrix offsets
        final int ido = dimTrait * iMatrix;
        final int jdo = dimTrait * jMatrix;

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

        final double[] diagQdi = new double[dimTrait];
        System.arraycopy(diagonalActualizations, ido, diagQdi, 0, dimTrait);
        final double[] diagQdj = new double[dimTrait];
        System.arraycopy(diagonalActualizations, jdo, diagQdj, 0, dimTrait);

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
//            final DenseMatrix64F QdiPip = matrix7;
//            CommonOps.multTransA(Qdi, Pip, QdiPip);

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
//            final DenseMatrix64F QdjPjp = matrix8;
//            CommonOps.multTransA(Qdj, Pjp, QdjPjp);

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
            diagonalDoubleProduct(Pip, diagQdi, QdiPipQdi);
            diagonalDoubleProduct(Pjp, diagQdj, QdjPjpQdj);
            CommonOps.add(QdiPipQdi, QdjPjpQdj, Pk);

            if (DEBUG) {
                System.err.println("Qdi: " + diagQdi);
//                System.err.println("\tQdiPip: " + QdiPip);
                System.err.println("\tQdiPipQdi: " + QdiPipQdi);
                System.err.println("\tQdj: " + diagQdj);
//                System.err.println("\tQdjPjp: " + QdjPjp);
                System.err.println("\tQdjPjpQdj: " + QdjPjpQdj);
            }

            // B. Partial mean

            if (TIMING) {
                endTime("peel3");
                startTime("peel4");
            }

            final double[] displacementi = vector1;
            final double[] displacementj = vector2;

//            final int ido = dimTrait * iMatrix;
//            final int jdo = dimTrait * jMatrix;

            for (int g = 0; g < dimTrait; ++g) {
                displacementi[g] = partials[ibo + g] - displacements[ido + g];
                displacementj[g] = partials[jbo + g] - displacements[jdo + g];
            }

//                final double[] tmp = new double[dimTrait];
            final double[] tmp = vector0;
            for (int g = 0; g < dimTrait; ++g) {
                double sum = 0.0;
                for (int h = 0; h < dimTrait; ++h) {
                    sum += diagQdi[g] * Pip.unsafe_get(g, h) * displacementi[h];
                    sum += diagQdj[g] * Pjp.unsafe_get(g, h) * displacementj[h];
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

    private static void  diagonalDoubleProduct(DenseMatrix64F source, double[] diagonalScales, DenseMatrix64F destination){
        double[] scales = new double[diagonalScales.length * diagonalScales.length];
        diagonalToMatrix(diagonalScales, scales);

        for (int i = 0; i < destination.data.length; ++i) {
                destination.data[i] = scales[i] * source.data[i];
        }
    }

    private static void  diagonalToMatrix(double[] diagonalScales, double[] scales) {
        int dim = diagonalScales.length;
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                scales[i * dim + j] = diagonalScales[i] * diagonalScales[j];
            }
        }
    }

    private double[] diagonalActualizations;
    protected double[] stationaryVariances;
//    private DenseMatrix64F matrix7;
//    private DenseMatrix64F matrix8;
    private DenseMatrix64F matrix9;
    private DenseMatrix64F matrix10;

    public double[] getStationaryVariance(int precisionIndex) {

        assert (stationaryVariances != null);

        final int offset = dimTrait * dimTrait * precisionIndex;

        double[] buffer = new double[dimTrait * dimTrait];

        System.arraycopy(stationaryVariances, offset, buffer, 0, dimTrait * dimTrait);

        return buffer;
    }
}

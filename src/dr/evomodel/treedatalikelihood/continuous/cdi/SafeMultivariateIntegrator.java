/*
 * SafeMultivariateIntegrator.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 *
 */

package dr.evomodel.treedatalikelihood.continuous.cdi;

import dr.math.MathUtils;
import dr.math.matrixAlgebra.WrappedVector;
import dr.math.matrixAlgebra.missingData.InversionResult;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

import static dr.math.matrixAlgebra.missingData.InversionResult.Code.*;
import static dr.math.matrixAlgebra.missingData.InversionResult.mult;
import static dr.math.matrixAlgebra.missingData.MissingOps.*;

/**
 * @author Marc A. Suchard
 */

public class SafeMultivariateIntegrator extends MultivariateIntegrator {

    private static final boolean DEBUG = false;
    private static final String ZERO_GLOBAL_REMAINDER_ADJOINT_PROPERTY =
            "beast.experimental.zeroGlobalRemainderAdjoint";

    public SafeMultivariateIntegrator(PrecisionType precisionType, int numTraits, int dimTrait, int dimProcess,
                                      int bufferCount, int diffusionCount) {
        super(precisionType, numTraits, dimTrait, dimProcess, bufferCount, diffusionCount);

        allocateStorage();

        effectiveDimensionOffset = PrecisionType.FULL.getEffectiveDimensionOffset(dimTrait);
        determinantOffset = PrecisionType.FULL.getDeterminantOffset(dimTrait);

        System.err.println("Trying SafeMultivariateIntegrator");
    }

    private void allocateStorage() {

        precisions = new double[dimTrait * dimTrait * bufferCount];
        variances = new double[dimTrait * dimTrait * bufferCount];

        vectorDelta = new double[dimTrait];
        vectorPMk = new double[dimTrait];
        matrixQjPjp = new DenseMatrix64F(dimTrait, dimTrait);
    }

    private static final boolean TIMING = false;

    @Override
    public void getBranchPrecision(int bufferIndex, int precisionIndex, DiffusionRepresentation precision) {
        getBranchPrecision(bufferIndex, ((DiffusionRepresentation.Dense) precision).matrix);
    }

    public void getBranchPrecision(int bufferIndex, double[] precision) {

        if (bufferIndex == -1) {
            throw new RuntimeException("Not yet implemented");
        }

        assert (precision != null);
        assert (precision.length >= dimTrait * dimTrait);

        System.arraycopy(precisions, bufferIndex * dimTrait * dimTrait,
                precision, 0, dimTrait * dimTrait);
    }

    @Override
    public void getBranchVariance(int bufferIndex, int precisionIndex, double[] precision) {
        getBranchVariance(bufferIndex, precision);
    }

    public void getBranchVariance(int bufferIndex, double[] variance) {

        if (bufferIndex == -1) {
            throw new RuntimeException("Not yet implemented");
        }

        assert (variance != null);
        assert (variance.length >= dimTrait * dimTrait);

        System.arraycopy(variances, bufferIndex * dimTrait * dimTrait,
                variance, 0, dimTrait * dimTrait);
    }

    @Override
    public void getRootPrecision(int priorBufferIndex, int precisionIndex, DiffusionRepresentation precision) {
        getRootPrecision(priorBufferIndex, ((DiffusionRepresentation.Dense) precision).matrix);
    }

    private void getRootPrecision(int priorBufferIndex, double[] precision) {

        assert (precision != null);
        assert (precision.length >= dimTrait * dimTrait);

        int priorOffset = dimPartial * priorBufferIndex;

        System.arraycopy(partials, priorOffset + dimTrait,
                precision, 0, dimTrait * dimTrait);
    }

    private double getEffectiveDimension(int iBuffer) {
        return partials[iBuffer * dimPartial + effectiveDimensionOffset];
    }

    @SuppressWarnings("unused")
    private void setEffectiveDimension(int iBuffer, double effDim) {
        partials[iBuffer * dimPartial + effectiveDimensionOffset] = effDim;
    }

    private double getPartialDeterminant(int iBuffer) {
        return partials[iBuffer * dimPartial + determinantOffset];
    }

    ///////////////////////////////////////////////////////////////////////////
    /// Setting variances, displacement and actualization vectors
    ///////////////////////////////////////////////////////////////////////////

    public void updateBrownianDiffusionMatrices(int precisionIndex, final int[] probabilityIndices,
                                                final double[] edgeLengths, final double[] driftRates,
                                                int updateCount) {

        super.updateBrownianDiffusionMatrices(precisionIndex, probabilityIndices, edgeLengths, driftRates, updateCount);

        assert (diffusions != null);
        assert (probabilityIndices.length >= updateCount);
        assert (edgeLengths.length >= updateCount);

        if (DEBUG) {
            System.err.println("Matrices (safe):");
        }

        final int matrixSize = dimProcess * dimProcess;
        final int unscaledOffset = matrixSize * precisionIndex;

        if (TIMING) {
            startTime("diffusion");
        }

        for (int up = 0; up < updateCount; ++up) {

            if (DEBUG) {
                System.err.println("\t" + probabilityIndices[up] + " <- " + edgeLengths[up]);
            }

            final double edgeLength = edgeLengths[up];

            final int scaledOffset = matrixSize * probabilityIndices[up];

            scale(diffusions, unscaledOffset, 1.0 / edgeLength, precisions, scaledOffset, matrixSize);
            scale(inverseDiffusions, unscaledOffset, edgeLength, variances, scaledOffset, matrixSize); // TODO Only if necessary
        }

        if (TIMING) {
            endTime("diffusion");
        }
    }

    static void scale(final double[] source,
                      final int sourceOffset,
                      final double scale,
                      final double[] destination,
                      final int destinationOffset,
                      final int length) {
        for (int i = 0; i < length; ++i) {
            destination[destinationOffset + i] = scale * source[sourceOffset + i];
        }
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

        // Read variance increments along descendant branches of k
        final DenseMatrix64F Vdi = wrap(variances, imo, dimTrait, dimTrait);
        final DenseMatrix64F Vdj = wrap(variances, jmo, dimTrait, dimTrait);

//        final DenseMatrix64F Pdi = wrap(precisions, imo, dimTrait, dimTrait); // TODO Only if needed
        final DenseMatrix64F Pdj = wrap(precisions, jmo, dimTrait, dimTrait); // TODO Only if needed

//        final DenseMatrix64F Vd = wrap(inverseDiffusions, precisionOffset, dimTrait, dimTrait);

        if (DEBUG) {
            System.err.println("updatePreOrderPartial for node " + iBuffer);
            System.err.println("\tVdj: " + Vdj);
            System.err.println("\tVdi: " + Vdi);
        }

        // For each trait // TODO in parallel
        for (int trait = 0; trait < numTraits; ++trait) {

            // A. Get current precision of k and j
            final DenseMatrix64F Pk = wrap(preOrderPartials, kbo + dimTrait, dimTrait, dimTrait);
//            final DenseMatrix64F Pj = wrap(partials, jbo + dimTrait, dimTrait, dimTrait);

//            final DenseMatrix64F Vk = wrap(preOrderPartials, kbo + dimTrait + dimTrait * dimTrait, dimTrait, dimTrait);
//            final DenseMatrix64F Vj = wrap(partials, jbo + dimTrait + dimTrait * dimTrait, dimTrait, dimTrait);

            // B. Inflate variance along sibling branch using matrix inversion
//            final DenseMatrix64F Vjp = matrix0;
            final DenseMatrix64F Pjp = matrixPjp;
            increaseVariances(jbo, jBuffer, Vdj, Pdj, Pjp, false);

            // Actualize
            final DenseMatrix64F QjPjp = matrixQjPjp;
            actualizePrecision(Pjp, QjPjp, jbo, jmo, jdo);

            // C. Compute prePartial mean
            final DenseMatrix64F Pip = matrixPip;
            CommonOps.add(Pk, Pjp, Pip);

            final DenseMatrix64F Vip = matrix1;
            safeInvertPrecision(Pip, Vip, false);

            final double[] delta = vectorDelta;
            computeDelta(jbo, jdo, delta);

//            final double[] tmp = vector0;
//            weightedAverage(preOrderPartials, kbo, Pk,
//                    delta, 0, QjPjp,
//                    preOrderPartials, ibo, Vip,
//                    dimTrait, tmp);
            safeWeightedAverage(
                    new WrappedVector.Raw(preOrderPartials, kbo, dimTrait),
                    Pk,
                    new WrappedVector.Raw(delta, 0, dimTrait),
                    QjPjp,
                    new WrappedVector.Raw(preOrderPartials, ibo, dimTrait),
                    Vip,
                    dimTrait);

            // Second pre-order channel: parent-side message for branch (v -> i), before branch transport.
            System.arraycopy(preOrderPartials, ibo, parentPreOrderPartials, ibo, dimTrait);
            unwrap(Pip, parentPreOrderPartials, ibo + dimTrait);
            unwrap(Vip, parentPreOrderPartials, ibo + dimTrait + dimTrait * dimTrait);

            scaleAndDriftMean(ibo, imo, ido);

            // C. Inflate variance along node branch
            final DenseMatrix64F Vi = Vip;
            actualizeVariance(Vip, ibo, imo, ido);
            inflateBranch(Vdi, Vip, Vi);

            final DenseMatrix64F Pi = matrixPk;
            safeInvert2(Vi, Pi, false);

            // X. Store precision results for node
            unwrap(Pi, preOrderPartials, ibo + dimTrait);
            unwrap(Vi, preOrderPartials, ibo + dimTrait + dimTrait * dimTrait);

            if (DEBUG) {
                System.err.println("trait: " + trait);
                System.err.println("pM: " + new WrappedVector.Raw(preOrderPartials, kbo, dimTrait));
                System.err.println("pP: " + Pk);
                System.err.println("sM: " + new WrappedVector.Raw(partials, jbo, dimTrait));
                DenseMatrix64F Pj = wrap(partials, jbo + dimTrait, dimTrait, dimTrait);
                DenseMatrix64F Vj = new DenseMatrix64F(dimTrait, dimTrait);
                CommonOps.invert(Pj, Vj);
                System.err.println("sP: " + Vj);
                System.err.println("sP: " + Pj);
                DenseMatrix64F Vjp = new DenseMatrix64F(dimTrait, dimTrait);
                CommonOps.invert(Pjp, Vjp);
                System.err.println("sVp: " + Vjp);
                System.err.println("sPp: " + Pjp);
                System.err.println("Pip: " + Pip);
                System.err.println("QiPip: " + QjPjp);
                System.err.println("cM: " + new WrappedVector.Raw(preOrderPartials, ibo, dimTrait));
                System.err.println("cV: " + Vi);
            }

            // Get ready for next trait
            kbo += dimPartialForTrait;
            ibo += dimPartialForTrait;
            jbo += dimPartialForTrait;
        }
    }

    private void inflateBranch(DenseMatrix64F Vj, DenseMatrix64F Vdj, DenseMatrix64F Vjp) {
        CommonOps.add(Vj, Vdj, Vjp);
    }

    void actualizePrecision(DenseMatrix64F P, DenseMatrix64F QP, int jbo, int jmo, int jdo) {
        CommonOps.scale(1.0, P, QP);
    }

    void actualizeVariance(DenseMatrix64F V, int ibo, int imo, int ido) {
        // Do nothing
    }

    void scaleAndDriftMean(int ibo, int imo, int ido) {
        // Do nothing
    }

    void computeDelta(int jbo, int jdo, double[] delta) {
        System.arraycopy(partials, jbo, delta, 0, dimTrait);
    }

    @Override
    protected void updatePartial(
            final int kBuffer,
            final int iBuffer,
            final int iMatrix,
            final int jBuffer,
            final int jMatrix,
            final boolean computeRemainders,
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

        // Read variance increments along descendant branches of k
        final DenseMatrix64F Vdi = wrap(variances, imo, dimTrait, dimTrait);
        final DenseMatrix64F Vdj = wrap(variances, jmo, dimTrait, dimTrait);

        final DenseMatrix64F Pdi = wrap(precisions, imo, dimTrait, dimTrait); // TODO Only if needed
        final DenseMatrix64F Pdj = wrap(precisions, jmo, dimTrait, dimTrait); // TODO Only if needed

        if (DEBUG) {
            System.err.println("variance diffusion: " + wrap(inverseDiffusions, precisionOffset, dimProcess, dimProcess));
            System.err.println("precisionOffset = " + precisionOffset);
            System.err.println("\tVdi: " + Vdi);
            System.err.println("\tVdj: " + Vdj);
        }

        // For each trait // TODO in parallel
        for (int trait = 0; trait < numTraits; ++trait) {

            // Layout, offset, dim
            // trait, 0, dT
            // precision, dT, dT * dT
            // variance, dT + dT * dT, dT * dT
            // scalar, dT + 2 * dT * dT, 1

            // Increase variance along the branches i -> k and j -> k

            final DenseMatrix64F Pip = matrixPip;
            final DenseMatrix64F Pjp = matrixPjp;


            InversionResult ci = increaseVariances(ibo, iBuffer, Vdi, Pdi, Pip, computeRemainders);
            InversionResult cj = increaseVariances(jbo, jBuffer, Vdj, Pdj, Pjp, computeRemainders);

            if (TIMING) {
                endTime("peel2");
                startTime("peel3");
            }

            // Compute partial mean and precision at node k

            // A. Partial precision and variance (for later use) using one matrix inversion
            final DenseMatrix64F Pk = matrixPk;
            computePartialPrecision(ido, jdo, imo, jmo, Pip, Pjp, Pk);

            if (TIMING) {
                endTime("peel3");
            }

            // B. Partial mean
            partialMean(ibo, jbo, kbo, ido, jdo);

            if (TIMING) {
                startTime("peel5");
            }

            // C. Store precision
            unwrap(Pk, partials, kbo + dimTrait);

            if (TIMING) {
                endTime("peel5");
            }

            if (DEBUG) {
                final DenseMatrix64F Pi = wrap(partials, ibo + dimTrait, dimTrait, dimTrait);
                final DenseMatrix64F Pj = wrap(partials, jbo + dimTrait, dimTrait, dimTrait);
                reportMeansAndPrecisions(trait, ibo, jbo, kbo, Pi, Pj, Pk);
            }

            // Computer remainder at node k
            double remainder = 0.0;

            if (computeRemainders) {

                if (DEBUG) {
                    reportInversions(ci, cj, Pip, Pjp);
                }

                if (TIMING) {
                    startTime("remain");
                }

                if (!(ci.getReturnCode() == NOT_OBSERVED || cj.getReturnCode() == NOT_OBSERVED)) {

                    // Inner products
                    double SS = computeSS(ibo, Pip, jbo, Pjp, kbo, Pk, dimTrait);

                    remainder += -0.5 * SS;

                    if (DEBUG) {
                        System.err.println("\t\t\tSS = " + (SS));
                    }
                } // End if remainder

                double effectiveDimension = getEffectiveDimension(iBuffer) + getEffectiveDimension(jBuffer);
                remainder += -effectiveDimension * LOG_SQRT_2_PI;

                double deti = 0;
                double detj = 0;
                if (!(ci.getReturnCode() == NOT_OBSERVED)) {
                    deti = ci.getLogDeterminant(); // TODO: for OU, use det(exp(M)) = exp(tr(M)) ? (Qdi = exp(-A l_i))
                }
                if (!(cj.getReturnCode() == NOT_OBSERVED)) {
                    detj = cj.getLogDeterminant();
                }
                remainder += -0.5 * (deti + detj);

                if (DEBUG) {
                    System.err.println("\t\t\tdeti = " + ci.getLogDeterminant());
                    System.err.println("\t\t\tdetj = " + cj.getLogDeterminant());
                    System.err.println("\t\tremainder: " + remainder);
                }

                if (TIMING) {
                    endTime("remain");
                }
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

    private void reportInversions(InversionResult ci, InversionResult cj,
                                  DenseMatrix64F Pip, DenseMatrix64F Pjp) {
        System.err.println("i status: " + ci);
        System.err.println("j status: " + cj);
        System.err.println("Pip: " + Pip);
        System.err.println("Pjp: " + Pjp);
    }

    private InversionResult increaseVariances(int ibo,
                                              int iBuffer,
                                              final DenseMatrix64F Vdi,
                                              final DenseMatrix64F Pdi,
                                              final DenseMatrix64F Pip,
                                              final boolean getDeterminant) {

        if (TIMING) {
            startTime("peel1");
        }

        // A. Get current precision of i and j
        final DenseMatrix64F Pi = wrap(partials, ibo + dimTrait, dimTrait, dimTrait);

        if (TIMING) {
            endTime("peel1");
            startTime("peel2");
        }

        // B. Integrate along branch using two matrix inversions

        final boolean useVariancei = anyDiagonalInfinities(Pi);
        InversionResult ci = null;

        if (useVariancei) {

            final DenseMatrix64F Vip = matrix0;
            final DenseMatrix64F Vi = wrap(partials, ibo + dimTrait + dimTrait * dimTrait, dimTrait, dimTrait);
            CommonOps.add(Vi, Vdi, Vip);
            if (allZeroOrInfinite(Vip)) {
                throw new RuntimeException("Zero-length branch on data is not allowed.");
            }
            ci = safeInvert2(Vip, Pip, getDeterminant);

        } else {

            final DenseMatrix64F tmp1 = matrix0;
            CommonOps.add(Pi, Pdi, tmp1);
            final DenseMatrix64F tmp2 = matrix1;
            safeInvertPrecision(tmp1, tmp2, false);
            CommonOps.mult(tmp2, Pi, tmp1);
            idMinusA(tmp1);
            if (getDeterminant) ci = safeDeterminant(tmp1, true);
            CommonOps.mult(Pi, tmp1, Pip);
            int effDim = (int) Math.round(getEffectiveDimension(iBuffer));
            if (getDeterminant && effDim > 0) {
                // effectiveDimension > 0 => a tip node & determinant not included earlier
                final InversionResult cP;
                double preCalculatedDeterminant = getPartialDeterminant(iBuffer);

                if (PrecisionType.FULL.isMissingDeterminantValue(preCalculatedDeterminant)) {
                    cP = safeDeterminant(Pi, true);
                } else {
                    InversionResult.Code code = InversionResult.getCode(dimTrait, effDim);
                    cP = new InversionResult(code, effDim, -preCalculatedDeterminant);
                }

                ci = mult(ci, cP);
            }
        }

        if (TIMING) {
            endTime("peel2");
        }

        return ci;
    }

    private static void idMinusA(DenseMatrix64F A) {
        CommonOps.scale(-1.0, A);
        for (int i = 0; i < A.numCols; i++) {
            A.set(i, i, 1.0 + A.get(i, i));
        }
    }

    private static boolean allZeroOrInfinite(DenseMatrix64F M) {
        for (int i = 0; i < M.getNumElements(); i++) {
            if (Double.isFinite(M.get(i)) && M.get(i) != 0.0) return false;
        }
        return true;
    }

    void computePartialPrecision(int ido, int jdo, int imo, int jmo,
                                 DenseMatrix64F Pip, DenseMatrix64F Pjp, DenseMatrix64F Pk) {
        CommonOps.add(Pip, Pjp, Pk);
    }

    void partialMean(int ibo, int jbo, int kbo,
                     int ido, int jdo) {
        if (TIMING) {
            startTime("peel4");
        }

        final double[] tmp = vectorPMk;
        weightedSum(partials, ibo, matrixPip, partials, jbo, matrixPjp, dimTrait, tmp);


        final WrappedVector kPartials = new WrappedVector.Raw(partials, kbo, dimTrait);
        final WrappedVector wrapTmp = new WrappedVector.Raw(tmp, 0, dimTrait);

        safeSolve(matrixPk, wrapTmp, kPartials, false);

        if (TIMING) {
            endTime("peel4");
            startTime("peel5");
        }
    }

    @Override
    public void calculateRootLogLikelihood(int rootBufferIndex, int priorBufferIndex, int precisionIndex,
                                           final double[] logLikelihoods,
                                           boolean incrementOuterProducts, boolean isIntegratedProcess) {
        assert (logLikelihoods.length == numTraits);

        assert (!incrementOuterProducts);

        updatePrecisionOffsetAndDeterminant(precisionIndex);

        if (DEBUG) {
            System.err.println("Root calculation for " + rootBufferIndex);
            System.err.println("Prior buffer index is " + priorBufferIndex);
        }

        int rootOffset = dimPartial * rootBufferIndex;
        int priorOffset = dimPartial * priorBufferIndex;

        // THIS IS THE INSTANTANEOUS PRECISION
        final DenseMatrix64F Pd = wrap(diffusions, precisionOffset, dimProcess, dimProcess); //TODO THIS IS THE INSTANTANEOUS PRECISION
//        final DenseMatrix64F Vd = wrap(inverseDiffusions, precisionOffset, dimTrait, dimTrait);

        // TODO For each trait in parallel
        for (int trait = 0; trait < numTraits; ++trait) {

            final DenseMatrix64F PPrior = wrap(partials, priorOffset + dimTrait, dimTrait, dimTrait);
            final DenseMatrix64F VPrior = wrap(partials, priorOffset + dimTrait + dimTrait * dimTrait, dimTrait, dimTrait);
// TODO VPRIOR IS JUST THE INVERSE OF PPRIOR
            if (DEBUG) {
//                System.out.println("SS root trait " + trait + ": " + Arrays.toString(PPrior.getData()));

            }

            // TODO Block below is for the conjugate prior ONLY
            {

                if (!isIntegratedProcess) {
                    final DenseMatrix64F PTmp = new DenseMatrix64F(dimTrait, dimTrait);
                    CommonOps.mult(Pd, PPrior, PTmp);
//                    CommonOps.mult(Pd, PPrior, PTmp);
//                    CommonOps.multTransB(PTmp, Pd, PPrior); //weighting by the instantaneous precision
                    PPrior.set(PTmp); // TODO What does this do?
                } else {
                    DenseMatrix64F Pdbis = new DenseMatrix64F(dimTrait, dimTrait);
                    blockUnwrap(Pd, Pdbis.data, 0, 0, 0, dimTrait);
                    blockUnwrap(Pd, Pdbis.data, dimProcess, dimProcess, 0, dimTrait);

                    final DenseMatrix64F PTmp = new DenseMatrix64F(dimTrait, dimTrait);
                    CommonOps.mult(Pdbis, PPrior, PTmp);
                    PPrior.set(PTmp);
                }
            }

            final DenseMatrix64F VTotal = new DenseMatrix64F(dimTrait, dimTrait);

            final DenseMatrix64F PTotal = new DenseMatrix64F(dimTrait, dimTrait);
            CommonOps.invert(VTotal, PTotal);  // TODO Does this do anything?
//            System.out.println("Vtotal " + trait + ": " + Arrays.toString(VTotal.getData()));
//            System.out.println("Ptotal " + trait + ": " + Arrays.toString(PTotal.getData()));
            InversionResult ctot = increaseVariances(rootOffset, rootBufferIndex, VPrior, PPrior, PTotal, true);

            double SS = weightedInnerProductOfDifferences(
                    partials, rootOffset,
                    partials, priorOffset,
                    PTotal, dimTrait);
            if (DEBUG) {
                System.out.println("SS root trait " + trait + ": " + SS);
            }

            double dettot = (ctot.getReturnCode() == NOT_OBSERVED) ? 0 : ctot.getLogDeterminant();
//            System.out.println("det root trait " + trait + ": " + dettot);
            if (DEBUG) {
//                System.err.println("partials: " + Arrays.toString(partials));
//                throw new RuntimeException("Stop here");
            }
            final double logLike = -0.5 * dettot - 0.5 * SS;

            final double remainder = remainders[rootBufferIndex * numTraits + trait];
//            System.out.println("remainder " + remainder);
            logLikelihoods[trait] = logLike + remainder;

            if (DEBUG) {
                System.err.print("mean:");
                for (int g = 0; g < dimTrait; ++g) {
                    System.err.print(" " + partials[rootOffset + g]);
                }
                System.err.println("");
                System.err.println("PRoot: " + wrap(partials, rootOffset + dimTrait, dimTrait, dimTrait));
                System.err.println("PPrior: " + PPrior);
                System.err.println("PTotal: " + PTotal);
                System.err.println("\n SS:" + SS);
                System.err.println("det:" + dettot);
                System.err.println("remainder:" + remainder);
                System.err.println("likelihood" + (logLike + remainder));
            }

            rootOffset += dimPartialForTrait;
            priorOffset += dimPartialForTrait;
        }

        if (DEBUG) {
            System.err.println("End");
        }
    }

//    private InversionResult computeBranchAdjustedPrecision(final double[] partials,
//                                                           final int bo,
//                                                           final DenseMatrix64F P,
//                                                           final DenseMatrix64F Pd,
//                                                           final DenseMatrix64F Vd,
//                                                           final double v,
//                                                           final DenseMatrix64F Pp) {
//        InversionResult c;
//        if (anyDiagonalInfinities(P)) {
//            // Inflate variance
//            final DenseMatrix64F Vp = matrix0;
//            final DenseMatrix64F Vi = wrap(partials, bo + dimTrait + dimTrait * dimTrait, dimTrait, dimTrait);
//
//            CommonOps.add(Vi, v, Vd, Vp);
//            c = safeInvert(Vp, Pp, true);
//        } else {
//            // Deflate precision
//            final DenseMatrix64F PPlusPd = matrix0;
//            CommonOps.add(P, 1.0 / v, Pd, PPlusPd);
//
//            final DenseMatrix64F PPlusPdInv = new DenseMatrix64F(dimTrait, dimTrait);
//            safeInvert(PPlusPd, PPlusPdInv, false);
//
//            CommonOps.mult(PPlusPdInv, P, Pp);
//            CommonOps.mult(P, Pp, PPlusPdInv);
//            CommonOps.add(P, -1, PPlusPdInv, Pp);
//            c = safeDeterminant(Pp, false);
//        }
//
//        return c;
//    }

    double computeSS(final int ibo,
                     final DenseMatrix64F Pip,
                     final int jbo,
                     final DenseMatrix64F Pjp,
                     final int kbo,
                     final DenseMatrix64F Pk,
                     final int dimTrait) {
        return weightedThreeInnerProductNormalized(partials, ibo, Pip,
                partials, jbo, Pjp,
                partials, kbo,
                vectorPMk, 0,
                dimTrait);
    }

    /**
     * Returns the displacement vector used in the remainder SS term for child i.
     * Base implementation uses the raw child mean; drift-aware subclasses override.
     */
    protected double[] getSSDisplacement(final int ibo, final int ido) {
        final double[] disp = new double[dimTrait];
        System.arraycopy(partials, ibo, disp, 0, dimTrait);
        return disp;
    }

    /**
     * Backward pass for one branch: computes bar(Vd_i), the adjoint of the branch covariance
     * entering {@link #increaseVariances(int, int, DenseMatrix64F, DenseMatrix64F, DenseMatrix64F, boolean)}.
     */
    @Override
    public void computeGlobalRemainderAdjointWrtBranchVariance(
            final int iBuffer,
            final int iMatrix,
            final int kBuffer,
            final double[] barVdiOut) {

        if (Boolean.getBoolean(ZERO_GLOBAL_REMAINDER_ADJOINT_PROPERTY)) {
            for (int i = 0; i < dimTrait * dimTrait; ++i) {
                barVdiOut[i] = 0.0;
            }
            return;
        }

        final double barSS = -0.5;
        final double barEll = -0.5;

        final int matrixOffset = dimTrait * dimTrait * iMatrix;
        final int displacementOffset = dimTrait * iMatrix;

        final DenseMatrix64F Vdi = new DenseMatrix64F(dimTrait, dimTrait);
        final DenseMatrix64F Pdi = new DenseMatrix64F(dimTrait, dimTrait);
        System.arraycopy(variances, matrixOffset, Vdi.data, 0, dimTrait * dimTrait);
        System.arraycopy(precisions, matrixOffset, Pdi.data, 0, dimTrait * dimTrait);

        final DenseMatrix64F barVdiAccum = new DenseMatrix64F(dimTrait, dimTrait);

        int ibo = dimPartial * iBuffer;
        int kbo = dimPartial * kBuffer;

        for (int trait = 0; trait < numTraits; ++trait) {
            final DenseMatrix64F Pip = new DenseMatrix64F(dimTrait, dimTrait);
            final InversionResult ci = increaseVariances(ibo, iBuffer, Vdi, Pdi, Pip, true);

            if (ci.getReturnCode() == NOT_OBSERVED) {
                ibo += dimPartialForTrait;
                kbo += dimPartialForTrait;
                continue;
            }

            final DenseMatrix64F barPip = computeBarPip(ibo, displacementOffset, kbo, barSS);
            accumulateIncreaseVariancesAdjoint(ibo, Pdi, Pip, barPip, barEll, barVdiAccum);

            ibo += dimPartialForTrait;
            kbo += dimPartialForTrait;
        }

        System.arraycopy(barVdiAccum.data, 0, barVdiOut, 0, dimTrait * dimTrait);
    }

    /**
     * bar(Pip) contribution from the SS term.
     */
    protected DenseMatrix64F computeBarPip(final int ibo,
                                           final int ido,
                                           final int kbo,
                                           final double barSS) {
        final double[] dispi = getSSDisplacement(ibo, ido);
        final double[] mk = new double[dimTrait];
        System.arraycopy(partials, kbo, mk, 0, dimTrait);

        final DenseMatrix64F barPip = new DenseMatrix64F(dimTrait, dimTrait);
        for (int g = 0; g < dimTrait; ++g) {
            final double rg = dispi[g] - mk[g];
            for (int h = 0; h < dimTrait; ++h) {
                barPip.set(g, h, barSS * rg * (dispi[h] - mk[h]));
            }
        }
        return barPip;
    }

    private void accumulateIncreaseVariancesAdjoint(final int ibo,
                                                    final DenseMatrix64F Pdi,
                                                    final DenseMatrix64F Pip,
                                                    final DenseMatrix64F barPip,
                                                    final double barEll,
                                                    final DenseMatrix64F barVdiAccum) {
        final DenseMatrix64F Pi = wrap(partials, ibo + dimTrait, dimTrait, dimTrait);
        final boolean useVariancei = anyDiagonalInfinities(Pi);

        final DenseMatrix64F tmp = new DenseMatrix64F(dimTrait, dimTrait);
        final DenseMatrix64F tmp2 = new DenseMatrix64F(dimTrait, dimTrait);

        if (useVariancei) {
            CommonOps.mult(Pip, barPip, tmp);
            CommonOps.mult(tmp, Pip, tmp2);
            for (int idx = 0; idx < dimTrait * dimTrait; ++idx) {
                barVdiAccum.data[idx] += -tmp2.data[idx] + barEll * Pip.data[idx];
            }
            return;
        }

        final DenseMatrix64F A = new DenseMatrix64F(dimTrait, dimTrait);
        CommonOps.add(Pi, Pdi, A);

        final DenseMatrix64F S = new DenseMatrix64F(dimTrait, dimTrait);
        CommonOps.invert(A, S);

        final DenseMatrix64F T = new DenseMatrix64F(dimTrait, dimTrait);
        CommonOps.mult(S, Pi, T);

        final DenseMatrix64F U = new DenseMatrix64F(dimTrait, dimTrait);
        for (int idx = 0; idx < dimTrait * dimTrait; ++idx) {
            U.data[idx] = -T.data[idx];
        }
        for (int g = 0; g < dimTrait; ++g) {
            U.add(g, g, 1.0);
        }

        final DenseMatrix64F UinvT = new DenseMatrix64F(dimTrait, dimTrait);
        CommonOps.invert(U, UinvT);
        CommonOps.transpose(UinvT);

        final DenseMatrix64F barU = new DenseMatrix64F(dimTrait, dimTrait);
        CommonOps.mult(Pi, barPip, barU);
        for (int idx = 0; idx < dimTrait * dimTrait; ++idx) {
            barU.data[idx] += barEll * UinvT.data[idx];
        }

        CommonOps.mult(S, barU, tmp);
        CommonOps.mult(tmp, Pi, tmp2);
        final DenseMatrix64F barPdi = new DenseMatrix64F(dimTrait, dimTrait);
        CommonOps.mult(tmp2, S, barPdi);

        CommonOps.mult(Pdi, barPdi, tmp);
        CommonOps.mult(tmp, Pdi, tmp2);
        for (int idx = 0; idx < dimTrait * dimTrait; ++idx) {
            barVdiAccum.data[idx] -= tmp2.data[idx];
        }
    }

    private final int effectiveDimensionOffset;
    private final int determinantOffset;

    private DenseMatrix64F matrixQjPjp;
    private double[] vectorDelta;
    double[] vectorPMk;
}

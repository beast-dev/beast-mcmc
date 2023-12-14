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
    public void getBranchPrecision(int bufferIndex, int precisionIndex, double[] precision) {
        getBranchPrecision(bufferIndex, precision);
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
    public void getRootPrecision(int priorBufferIndex, int precisionIndex, double[] precision) {
        getRootPrecision(priorBufferIndex, precision);
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

        final DenseMatrix64F Pd = wrap(diffusions, precisionOffset, dimProcess, dimProcess);
//        final DenseMatrix64F Vd = wrap(inverseDiffusions, precisionOffset, dimTrait, dimTrait);

        // TODO For each trait in parallel
        for (int trait = 0; trait < numTraits; ++trait) {

            final DenseMatrix64F PPrior = wrap(partials, priorOffset + dimTrait, dimTrait, dimTrait);
            final DenseMatrix64F VPrior = wrap(partials, priorOffset + dimTrait + dimTrait * dimTrait, dimTrait, dimTrait);


            // TODO Block below is for the conjugate prior ONLY
            {

                if (!isIntegratedProcess) {
                    final DenseMatrix64F PTmp = new DenseMatrix64F(dimTrait, dimTrait);
                    CommonOps.mult(Pd, PPrior, PTmp);
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

            InversionResult ctot = increaseVariances(rootOffset, rootBufferIndex, VPrior, PPrior, PTotal, true);

            double SS = weightedInnerProductOfDifferences(
                    partials, rootOffset,
                    partials, priorOffset,
                    PTotal, dimTrait);

            double dettot = (ctot.getReturnCode() == NOT_OBSERVED) ? 0 : ctot.getLogDeterminant();

            final double logLike = -0.5 * dettot - 0.5 * SS;

            final double remainder = remainders[rootBufferIndex * numTraits + trait];
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

    private final int effectiveDimensionOffset;
    private final int determinantOffset;

    private DenseMatrix64F matrixQjPjp;
    private double[] vectorDelta;
    double[] vectorPMk;
}

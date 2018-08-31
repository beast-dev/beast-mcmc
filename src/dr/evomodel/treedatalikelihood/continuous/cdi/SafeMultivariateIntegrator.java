package dr.evomodel.treedatalikelihood.continuous.cdi;

import dr.evomodel.treedatalikelihood.PostOrderStatistics;
import dr.math.matrixAlgebra.WrappedVector;
import dr.math.matrixAlgebra.missingData.InversionResult;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

import static dr.math.matrixAlgebra.missingData.InversionResult.Code.NOT_OBSERVED;
import static dr.math.matrixAlgebra.missingData.MissingOps.*;

/**
 * @author Marc A. Suchard
 */

public class SafeMultivariateIntegrator extends MultivariateIntegrator {

    private static boolean DEBUG = false;

    public SafeMultivariateIntegrator(PrecisionType precisionType, int numTraits, int dimTrait, int bufferCount,
                                      int diffusionCount) {
        super(precisionType, numTraits, dimTrait, bufferCount, diffusionCount);

        allocateStorage();

        System.err.println("Trying SafeMultivariateIntegrator");
    }

    private void allocateStorage() {

        precisions = new double[dimTrait * dimTrait * bufferCount];
        variances = new double[dimTrait * dimTrait * bufferCount];

        vectorDelta = new double[dimTrait];

        matrixQjPjp = new DenseMatrix64F(dimTrait, dimTrait);
    }

    private static final boolean TIMING = false;

    @Override
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
    public void getRootPrecision(int priorBufferIndex, double[] precision) {

        assert (precision != null);
        assert (precision.length >= dimTrait * dimTrait);

        int priorOffset = dimPartial * priorBufferIndex;

        System.arraycopy(partials, priorOffset + dimTrait,
                precision, 0, dimTrait * dimTrait);
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

        final int matrixSize = dimTrait * dimTrait;
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
            final DenseMatrix64F Vjp = matrix0;
            final DenseMatrix64F Pjp = matrixPjp;
            increaseVariances(jbo, Vdj, Pdj, Pjp, false);

            // Actualize
            final DenseMatrix64F QjPjp = matrixQjPjp;
            actualizePrecision(Pjp, QjPjp, jbo, jmo, jdo);

            // C. Compute prePartial mean
            final DenseMatrix64F Pip = matrixPip;
            CommonOps.add(Pk, Pjp, Pip);

            final DenseMatrix64F Vip = matrix1;
            safeInvert(Pip, Vip, false);

            final double[] delta = vectorDelta;
            computeDelta(jbo, jdo, delta);

            final double[] tmp = vector0;
            weightedAverage(preOrderPartials, kbo, Pk,
                    delta, 0, QjPjp,
                    preOrderPartials, ibo, Vip,
                    dimTrait, tmp);

            scaleAndDriftMean(ibo, imo, ido);

            // C. Inflate variance along node branch
            final DenseMatrix64F Vi = Vip;
            actualizeVariance(Vip, ibo, imo, ido);
            inflateBranch(Vdi, Vip, Vi);

            final DenseMatrix64F Pi = matrixPk;
            safeInvert(Vi, Pi, false);

            // X. Store precision results for node
            unwrap(Pi, preOrderPartials, ibo + dimTrait);
            unwrap(Vi, preOrderPartials, ibo + dimTrait + dimTrait * dimTrait);

            if (DEBUG) {
                System.err.println("trait: " + trait);
                System.err.println("pM: " + new WrappedVector.Raw(preOrderPartials, kbo, dimTrait));
                System.err.println("pP: " + Pk);
                System.err.println("sM: " + new WrappedVector.Raw(partials, jbo, dimTrait));
                DenseMatrix64F Pj = wrap(partials, ibo + dimTrait, dimTrait, dimTrait);
                DenseMatrix64F Vj = new DenseMatrix64F(dimTrait, dimTrait);
                CommonOps.invert(Pj, Vj);
                System.err.println("sP: " + Vj);
                System.err.println("sP: " + Pj);
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
            final PostOrderStatistics.Continuous statistics) {

        if (statistics == PostOrderStatistics.Continuous.REMAINDERS_AND_WISHART) {
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

        final DenseMatrix64F Vd = wrap(inverseDiffusions, precisionOffset, dimTrait, dimTrait);

        if (DEBUG) {
            System.err.println("variance diffusion: " + Vd);
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


            InversionResult ci = increaseVariances(ibo, Vdi, Pdi, Pip, true);
            InversionResult cj = increaseVariances(jbo, Vdj, Pdj, Pjp, true);

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
            InversionResult ck = partialMean(ibo, jbo, kbo, ido, jdo);

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

            if (DEBUG) {
                reportInversions(ci, cj, ck, Pip, Pjp);
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

            int dimensionChange = ci.getEffectiveDimension() + cj.getEffectiveDimension()
                    - ck.getEffectiveDimension();

            remainder += -dimensionChange * LOG_SQRT_2_PI;

            double deti = 0;
            double detj = 0;
            double detk = 0;
            if (!(ci.getReturnCode() == NOT_OBSERVED)) {
                deti = Math.log(ci.getDeterminant()); // TODO: for OU, use det(exp(M)) = exp(tr(M)) ? (Qdi = exp(-A l_i))
            }
            if (!(cj.getReturnCode() == NOT_OBSERVED)) {
                detj = Math.log(cj.getDeterminant());
            }
            if (!(ck.getReturnCode() == NOT_OBSERVED)) {
                detk = Math.log(ck.getDeterminant());
            }
            remainder += -0.5 * (deti + detj + detk);

            // TODO Can get SSi + SSj - SSk from inner product w.r.t Pt (see outer-products below)?

            if (DEBUG) {
                System.err.println("\t\t\tdeti = " + Math.log(ci.getDeterminant()));
                System.err.println("\t\t\tdetj = " + Math.log(cj.getDeterminant()));
                System.err.println("\t\t\tdetk = " + Math.log(ck.getDeterminant()));
                System.err.println("\t\tremainder: " + remainder);
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

    private void reportInversions(InversionResult ci, InversionResult cj, InversionResult ck,
                                  DenseMatrix64F Pip, DenseMatrix64F Pjp) {
        System.err.println("i status: " + ci);
        System.err.println("j status: " + cj);
        System.err.println("k status: " + ck);
        System.err.println("Pip: " + Pip);
        System.err.println("Pjp: " + Pjp);
    }

    private InversionResult increaseVariances(int ibo,
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
//                CommonOps.add(Vi, vi, Vd, Vip);  // TODO Fix
            CommonOps.add(Vi, Vdi, Vip);
            ci = safeInvert(Vip, Pip, getDeterminant);

        } else {

            final DenseMatrix64F PiPlusPd = matrix0;
//                CommonOps.add(Pi, 1.0 / vi, Pd, PiPlusPd); // TODO Fix
            CommonOps.add(Pi, Pdi, PiPlusPd);
            final DenseMatrix64F PiPlusPdInv = new DenseMatrix64F(dimTrait, dimTrait);
            safeInvert(PiPlusPd, PiPlusPdInv, false);
            CommonOps.mult(PiPlusPdInv, Pi, Pip);
            CommonOps.mult(Pi, Pip, PiPlusPdInv);
            CommonOps.add(Pi, -1, PiPlusPdInv, Pip);
            if (getDeterminant) ci = safeDeterminant(Pip, false);
        }

        if (TIMING) {
            endTime("peel2");
        }

        return ci;
    }

    void computePartialPrecision(int ido, int jdo, int imo, int jmo,
                                 DenseMatrix64F Pip, DenseMatrix64F Pjp, DenseMatrix64F Pk) {
        CommonOps.add(Pip, Pjp, Pk);
    }

    InversionResult partialMean(int ibo, int jbo, int kbo,
                                int ido, int jdo) {
        if (TIMING) {
            startTime("peel4");
        }

        final double[] tmp = vector0;
        weightedSum(partials, ibo, matrixPip, partials, jbo, matrixPjp, dimTrait, tmp);


        final WrappedVector kPartials = new WrappedVector.Raw(partials, kbo, dimTrait);
        final WrappedVector wrapTmp = new WrappedVector.Raw(tmp, 0, dimTrait);

        InversionResult ck = safeSolve(matrixPk, wrapTmp, kPartials, true);

        if (TIMING) {
            endTime("peel4");
            startTime("peel5");
        }
        return ck;
    }


//    private final Map<String, Long> startTimes = new HashMap<String, Long>();
//
//    private void startTime(String key) {
//        startTimes.put(key, System.nanoTime());
//    }
//
//    private void endTime(String key) {
//        long start = startTimes.get(key);
//
//        Long total = times.get(key);
//        if (total == null) {
//            total = new Long(0);
//        }
//
//        long run = total + (System.nanoTime() - start);
//        times.put(key, run);
//
////            System.err.println("run = " + run);
////            System.exit(-1);
//    }

//        private void incrementTiming(long start, long end, String key) {
//            Long total = times.get(key);
//
//            System.err.println(start + " " + end + " " + key);
//            System.exit(-1);
//            if (total == null) {
//                total = new Long(0);
//                times.put(key, total);
//            }
//            total += (end - start);
////            times.put(key, total);
//        }

    @Override
    public void calculateRootLogLikelihood(int rootBufferIndex, int priorBufferIndex, final double[] logLikelihoods,
                                           PostOrderStatistics.Continuous statistics) {

        assert (logLikelihoods.length == numTraits);

        assert (statistics != PostOrderStatistics.Continuous.REMAINDERS_AND_WISHART);

        if (DEBUG) {
            System.err.println("Root calculation for " + rootBufferIndex);
            System.err.println("Prior buffer index is " + priorBufferIndex);
        }

        int rootOffset = dimPartial * rootBufferIndex;
        int priorOffset = dimPartial * priorBufferIndex;

        final DenseMatrix64F Pd = wrap(diffusions, precisionOffset, dimTrait, dimTrait);
//        final DenseMatrix64F Vd = wrap(inverseDiffusions, precisionOffset, dimTrait, dimTrait);

        // TODO For each trait in parallel
        for (int trait = 0; trait < numTraits; ++trait) {

            final DenseMatrix64F PRoot = wrap(partials, rootOffset + dimTrait, dimTrait, dimTrait);
            final DenseMatrix64F PPrior = wrap(partials, priorOffset + dimTrait, dimTrait, dimTrait);

            // TODO Block below is for the conjugate prior ONLY
            {

                final DenseMatrix64F PTmp = new DenseMatrix64F(dimTrait, dimTrait);
                CommonOps.mult(Pd, PPrior, PTmp);
                PPrior.set(PTmp); // TODO What does this do?
            }

            final DenseMatrix64F VTotal = new DenseMatrix64F(dimTrait, dimTrait);
//            CommonOps.add(VRoot, VPrior, VTotal);

            final DenseMatrix64F PTotal = new DenseMatrix64F(dimTrait, dimTrait);
            CommonOps.invert(VTotal, PTotal);  // TODO Does this do anything?

            final boolean useVariance = anyDiagonalInfinities(PRoot);
            InversionResult ctot;

            if (useVariance) {
                final DenseMatrix64F Vroot = wrap(partials, rootOffset + dimTrait + dimTrait * dimTrait, dimTrait, dimTrait);
                final DenseMatrix64F VPrior = wrap(partials, priorOffset + dimTrait + dimTrait * dimTrait, dimTrait, dimTrait);
                CommonOps.add(Vroot, VPrior, VTotal);
                ctot = safeInvert(VTotal, PTotal, true);
            } else {
                final DenseMatrix64F tmp1 = new DenseMatrix64F(dimTrait, dimTrait);
                final DenseMatrix64F tmp2 = new DenseMatrix64F(dimTrait, dimTrait);
                CommonOps.add(PRoot, PPrior, PTotal);
                safeInvert(PTotal, VTotal, false);
                CommonOps.mult(VTotal, PRoot, tmp1);
                CommonOps.mult(PRoot, tmp1, tmp2);
                CommonOps.add(PRoot, -1.0, tmp2, PTotal);
                ctot = safeDeterminant(PTotal, false);
            }


//            double SS = 0;
//            for (int g = 0; g < dimTrait; ++g) {
//                final double gDifference = partials[rootOffset + g] - partials[priorOffset + g];
//
//                for (int h = 0; h < dimTrait; ++h) {
//                    final double hDifference = partials[rootOffset + h] - partials[priorOffset + h];
//
//                    SS += gDifference * PTotal.unsafe_get(g, h) * hDifference;
//                }
//            }
            double SS = weightedInnerProductOfDifferences(
                    partials, rootOffset,
                    partials, priorOffset,
                    PTotal, dimTrait);

            final double logLike = -ctot.getEffectiveDimension() * LOG_SQRT_2_PI
//                    - 0.5 * Math.log(CommonOps.det(VTotal))
//                    + 0.5 * Math.log(CommonOps.det(PTotal))
                    - 0.5 * Math.log(ctot.getDeterminant())
                    - 0.5 * SS;

            final double remainder = remainders[rootBufferIndex * numTraits + trait];
            logLikelihoods[trait] = logLike + remainder;

//            if (incrementOuterProducts) {
//
//                assert false : "Should not get here";
//
////                int opo = dimTrait * dimTrait * trait;
////                int opd = precisionOffset;
////
////                double rootScalar = partials[rootOffset + dimTrait + 2 * dimTrait * dimTrait];
////                final double priorScalar = partials[priorOffset + dimTrait];
////
////                if (!Double.isInfinite(priorScalar)) {
////                    rootScalar = rootScalar * priorScalar / (rootScalar + priorScalar);
////                }
////
////                for (int g = 0; g < dimTrait; ++g) {
////                    final double gDifference = partials[rootOffset + g] - partials[priorOffset + g];
////
////                    for (int h = 0; h < dimTrait; ++h) {
////                        final double hDifference = partials[rootOffset + h] - partials[priorOffset + h];
////
////                        outerProducts[opo] += gDifference * hDifference
//////                                    * PTotal.unsafe_get(g, h) / diffusions[opd];
////                                * rootScalar;
////                        ++opo;
////                        ++opd;
////                    }
////                }
////
////                degreesOfFreedom[trait] += 1; // increment degrees-of-freedom
//            }

            if (DEBUG) {
                System.err.print("mean:");
                for (int g = 0; g < dimTrait; ++g) {
                    System.err.print(" " + partials[rootOffset + g]);
                }
                System.err.println("");
                System.err.println("PRoot: " + PRoot);
                System.err.println("PPrior: " + PPrior);
                System.err.println("PTotal: " + PTotal);
                System.err.println("\t" + logLike + " " + (logLike + remainder));
                System.err.println("\t" + remainder);

//                if (incrementOuterProducts) {
//                    System.err.println("Outer-products:" + wrap(outerProducts, dimTrait * dimTrait * trait, dimTrait, dimTrait));
//                }
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
        return weightedThreeInnerProduct(partials, ibo, Pip,
                partials, jbo, Pjp,
                partials, kbo, Pk,
                dimTrait);
    }

    private DenseMatrix64F matrixQjPjp;
    private double[] vectorDelta;
}

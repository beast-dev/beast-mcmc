package dr.evomodel.treedatalikelihood.continuous.cdi;

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
        
        System.err.println("Trying SafeMultivariateIntegrator");
    }

    private static final boolean TIMING = false;

//    @Override
//    public void updatePreOrderPartial(
//            final int kBuffer, // parent
//            final int iBuffer, // node
//            final int iMatrix,
//            final int jBuffer, // sibling
//            final int jMatrix) {
//
//        // Determine buffer offsets
//        int kbo = dimPartial * kBuffer;
//        int ibo = dimPartial * iBuffer;
//        int jbo = dimPartial * jBuffer;
//
//        // Determine matrix offsets
//        final int imo = dimMatrix * iMatrix;
//        final int jmo = dimMatrix * jMatrix;
//
//        // Read variance increments along descendant branches of k
//        final double vi = branchLengths[imo];
//        final double vj = branchLengths[jmo];
//
//        final DenseMatrix64F Vd = wrap(inverseDiffusions, precisionOffset, dimTrait, dimTrait);
//
//        if (DEBUG) {
//            System.err.println("updatePreOrderPartial for node " + iBuffer);
////                System.err.println("variance diffusion: " + Vd);
//            System.err.println("\tvi: " + vi + " vj: " + vj);
////                System.err.println("precisionOffset = " + precisionOffset);
//        }
//
//        // For each trait // TODO in parallel
//        for (int trait = 0; trait < numTraits; ++trait) {
//
//            // A. Get current precision of k and j
//            final DenseMatrix64F Pk = wrap(prePartials, kbo + dimTrait, dimTrait, dimTrait);
////                final DenseMatrix64F Pj = wrap(partials, jbo + dimTrait, dimTrait, dimTrait);
//
////                final DenseMatrix64F Vk = wrap(prePartials, kbo + dimTrait + dimTrait * dimTrait, dimTrait, dimTrait);
//            final DenseMatrix64F Vj = wrap(partials, jbo + dimTrait + dimTrait * dimTrait, dimTrait, dimTrait);
//
//            // B. Inflate variance along sibling branch using matrix inversion
////                final DenseMatrix64F Vjp = new DenseMatrix64F(dimTrait, dimTrait);
//            final DenseMatrix64F Vjp = matrix0;
//            CommonOps.add(Vj, vj, Vd, Vjp);
//
////                final DenseMatrix64F Pjp = new DenseMatrix64F(dimTrait, dimTrait);
//            final DenseMatrix64F Pjp = matrix1;
//            InversionResult cj = safeInvert(Vjp, Pjp, false);
//
////                final DenseMatrix64F Pip = new DenseMatrix64F(dimTrait, dimTrait);
//            final DenseMatrix64F Pip = matrix2;
//            CommonOps.add(Pk, Pjp, Pip);
//
////                final DenseMatrix64F Vip = new DenseMatrix64F(dimTrait, dimTrait);
//            final DenseMatrix64F Vip = matrix3;
//            InversionResult cip = safeInvert(Pip, Vip, false);
//
//            // C. Compute prePartial mean
////                final double[] tmp = new double[dimTrait];
//            final double[] tmp = vector0;
//            for (int g = 0; g < dimTrait; ++g) {
//                double sum = 0.0;
//                for (int h = 0; h < dimTrait; ++h) {
//                    sum += Pk.unsafe_get(g, h) * prePartials[kbo + h]; // Read parent
//                    sum += Pjp.unsafe_get(g, h) * partials[jbo + h];   // Read sibling
//                }
//                tmp[g] = sum;
//            }
//            for (int g = 0; g < dimTrait; ++g) {
//                double sum = 0.0;
//                for (int h = 0; h < dimTrait; ++h) {
//                    sum += Vip.unsafe_get(g, h) * tmp[h];
//                }
//                prePartials[ibo + g] = sum; // Write node
//            }
//
//            // C. Inflate variance along node branch
//            final DenseMatrix64F Vi = Vip;
//            CommonOps.add(vi, Vd, Vip, Vi);
//
////                final DenseMatrix64F Pi = new DenseMatrix64F(dimTrait, dimTrait);
//            final DenseMatrix64F Pi = matrix4;
//            InversionResult ci = safeInvert(Vi, Pi, false);
//
//            // X. Store precision results for node
//            unwrap(Pi, prePartials, ibo + dimTrait);
//            unwrap(Vi, prePartials, ibo + dimTrait + dimTrait * dimTrait);
//
//            if (DEBUG) {
//                System.err.println("trait: " + trait);
//                System.err.println("pM: " + new WrappedVector.Raw(prePartials, kbo, dimTrait));
//                System.err.println("pP: " + Pk);
//                System.err.println("sM: " + new WrappedVector.Raw(partials, jbo, dimTrait));
//                System.err.println("sV: " + Vj);
//                System.err.println("sVp: " + Vjp);
//                System.err.println("sPp: " + Pjp);
//                System.err.println("Pip: " + Pip);
//                System.err.println("cM: " + new WrappedVector.Raw(prePartials, ibo, dimTrait));
//                System.err.println("cV: " + Vi);
//            }
//
//            // Get ready for next trait
//            kbo += dimPartialForTrait;
//            ibo += dimPartialForTrait;
//            jbo += dimPartialForTrait;
//        }
//    }

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
        final int imo = dimMatrix * iMatrix;
        final int jmo = dimMatrix * jMatrix;

        // Read variance increments along descendant branches of k
        final double vi = branchLengths[imo];
        final double vj = branchLengths[jmo];

        final DenseMatrix64F Vd = wrap(inverseDiffusions, precisionOffset, dimTrait, dimTrait);
        final DenseMatrix64F Pd = wrap(diffusions, precisionOffset, dimTrait, dimTrait);

        if (DEBUG) {
            System.err.println("variance diffusion: " + Vd);
            System.err.println("\tvi: " + vi + " vj: " + vj);
            System.err.println("precisionOffset = " + precisionOffset);
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
            final double lpi = partials[ibo + dimTrait + 2 * dimTrait * dimTrait];
            final double lpj = partials[jbo + dimTrait + 2 * dimTrait * dimTrait];

            final DenseMatrix64F Pi = wrap(partials, ibo + dimTrait, dimTrait, dimTrait);
            final DenseMatrix64F Pj = wrap(partials, jbo + dimTrait, dimTrait, dimTrait);

            if (TIMING) {
                endTime("peel1");
                startTime("peel2");
            }

            // B. Integrate along branch using two matrix inversions
            @SuppressWarnings("SpellCheckingInspection")
            final double lpip = Double.isInfinite(lpi) ?
                    1.0 / vi : lpi / (1.0 + lpi * vi);
            @SuppressWarnings("SpellCheckingInspection")
            final double lpjp = Double.isInfinite(lpj) ?
                    1.0 / vj : lpj / (1.0 + lpj * vj);

            final DenseMatrix64F Pip = matrix2;
            final DenseMatrix64F Pjp = matrix3;

            InversionResult ci;
            InversionResult cj;

            final boolean iUseVariance = anyDiagonalInfinities(Pi);
            final boolean jUseVariance = anyDiagonalInfinities(Pj);

            if (iUseVariance) {

                final DenseMatrix64F Vip = matrix0;
                final DenseMatrix64F Vi = wrap(partials, ibo + dimTrait + dimTrait * dimTrait, dimTrait, dimTrait);
                CommonOps.add(Vi, vi, Vd, Vip);
                ci = safeInvert(Vip, Pip, true);

            } else {

                final DenseMatrix64F PiPlusPd = matrix0;
                CommonOps.add(Pi, 1.0 / vi, Pd, PiPlusPd);
                final DenseMatrix64F PiPlusPdInv = new DenseMatrix64F(dimTrait, dimTrait);
                safeInvert(PiPlusPd, PiPlusPdInv, false);
                CommonOps.mult(PiPlusPdInv, Pi, Pip);
                CommonOps.mult(Pi, Pip, PiPlusPdInv);
                CommonOps.add(Pi, -1, PiPlusPdInv, Pip);
                ci = safeDeterminant(Pip, false);
            }
            // Or ...
//            InversionResult ci = computeBranchAdjustedPrecision(partials, ibo, Pi, Pd, Vd, vi, Pip);

            if (jUseVariance) {

                final DenseMatrix64F Vjp = matrix1;
                final DenseMatrix64F Vj = wrap(partials, jbo + dimTrait + dimTrait * dimTrait, dimTrait, dimTrait);
                CommonOps.add(Vj, vj, Vd, Vjp);
                cj = safeInvert(Vjp, Pjp, true);

            } else {

                final DenseMatrix64F PjPlusPd = matrix1;
                CommonOps.add(Pj, 1.0 / vj, Pd, PjPlusPd);
                final DenseMatrix64F PjPlusPdInv = new DenseMatrix64F(dimTrait, dimTrait);
                safeInvert(PjPlusPd, PjPlusPdInv, false);
                CommonOps.mult(PjPlusPdInv, Pj, Pjp);
                CommonOps.mult(Pj, Pjp, PjPlusPdInv);
                CommonOps.add(Pj, -1, PjPlusPdInv, Pjp);
                cj = safeDeterminant(Pjp, false);
            }
            // Or ...
//            InversionResult cj = computeBranchAdjustedPrecision(partials,jbo, Pj, Pd, Vd, vj, Pjp);

            if (TIMING) {
                endTime("peel2");
                startTime("peel2a");
            }

            if (TIMING) {
                endTime("peel2a");
                startTime("peel3");
            }

            // Compute partial mean and precision at node k

            // A. Partial precision and variance (for later use) using one matrix inversion
            final double lpk = lpip + lpjp;

            final DenseMatrix64F Pk = matrix4;
            CommonOps.add(Pip, Pjp, Pk);

            // B. Partial mean

            if (TIMING) {
                endTime("peel3");
                startTime("peel4");
            }

//                final double[] tmp = new double[dimTrait];
            final double[] tmp = vector0;
//            for (int g = 0; g < dimTrait; ++g) {
//                double sum = 0.0;
//                for (int h = 0; h < dimTrait; ++h) {
//                    sum += Pip.unsafe_get(g, h) * partials[ibo + h];
//                    sum += Pjp.unsafe_get(g, h) * partials[jbo + h];
//                }
//                tmp[g] = sum;
//            }
            weightedSum(partials, ibo, Pip, partials, jbo, Pjp, dimTrait, tmp);


//            for (int g = 0; g < dimTrait; ++g) {
//                double sum = 0.0;
//                for (int h = 0; h < dimTrait; ++h) {
//                    sum += Vk.unsafe_get(g, h) * tmp[h];
//                }
//                partials[kbo + g] = sum;
//            }

            final WrappedVector kPartials = new WrappedVector.Raw(partials, kbo, dimTrait);
            final WrappedVector wrapTmp = new WrappedVector.Raw(tmp, 0, dimTrait);

            InversionResult ck = safeSolve(Pk, wrapTmp, kPartials, true);

            if (TIMING) {
                endTime("peel4");
                startTime("peel5");
            }

            // C. Store precision
            partials[kbo + dimTrait + 2 * dimTrait * dimTrait] = lpk;

            unwrap(Pk, partials, kbo + dimTrait);

            if (TIMING) {
                endTime("peel5");
            }

            if (DEBUG) {
                reportMeansAndPrecisions(trait, ibo, jbo, kbo, Pi, Pj, Pk);
            }

            // Computer remainder at node k
            double remainder = 0.0;

            if (DEBUG) {
                reportInversions(ci, cj, ck, Pip, Pjp);

            }

            if (!(ci.getReturnCode() == NOT_OBSERVED || cj.getReturnCode() == NOT_OBSERVED)) {

                if (TIMING) {
                    startTime("remain");
                }

                // Inner products
                double SS = weightedThreeInnerProduct(
                        partials, ibo, Pip,
                        partials, jbo, Pjp,
                        partials, kbo, Pk,
                        dimTrait);

                int dimensionChange = ci.getEffectiveDimension() + cj.getEffectiveDimension()
                        - ck.getEffectiveDimension();

                remainder += -dimensionChange * LOG_SQRT_2_PI - 0.5 *
                        (Math.log(ci.getDeterminant()) + Math.log(cj.getDeterminant()) + Math.log(ck.getDeterminant()))
                        - 0.5 * SS;

                // TODO Can get SSi + SSj - SSk from inner product w.r.t Pt (see outer-products below)?

                if (DEBUG) {
                    System.err.println("\t\t\tDetI = " + Math.log(ci.getDeterminant()));
                    System.err.println("\t\t\tDetJ = " + Math.log(ci.getDeterminant()));
                    System.err.println("\t\t\tDetK = " + Math.log(ci.getDeterminant()));
                    System.err.println("\t\tremainder: " + remainder);
                }

                if (TIMING) {
                    endTime("remain");
                }

            } // End if remainder

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

    void reportInversions(InversionResult ci, InversionResult cj, InversionResult ck,
                          DenseMatrix64F Pip, DenseMatrix64F Pjp) {
        System.err.println("i status: " + ci);
        System.err.println("j status: " + cj);
        System.err.println("k status: " + ck);
        System.err.println("Pip: " + Pip);
        System.err.println("Pjp: " + Pjp);
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
                                           boolean incrementOuterProducts) {
        assert(logLikelihoods.length == numTraits);

        assert(!incrementOuterProducts);

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

            final DenseMatrix64F tmp1 = new DenseMatrix64F(dimTrait, dimTrait);
            final DenseMatrix64F tmp2 = new DenseMatrix64F(dimTrait, dimTrait);
            CommonOps.add(PRoot, PPrior, PTotal);
            CommonOps.invert(PTotal, VTotal);
            CommonOps.mult(VTotal, PRoot, tmp1);
            CommonOps.mult(PRoot, tmp1, tmp2);
            CommonOps.add(PRoot, -1.0, tmp2, PTotal);

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

            final double logLike = -dimTrait * LOG_SQRT_2_PI
//                    - 0.5 * Math.log(CommonOps.det(VTotal))
                    + 0.5 * Math.log(CommonOps.det(PTotal))
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

    private InversionResult computeBranchAdjustedPrecision(final double[] partials,
                                                           final int bo,
                                                           final DenseMatrix64F P,
                                                           final DenseMatrix64F Pd,
                                                           final DenseMatrix64F Vd,
                                                           final double v,
                                                           final DenseMatrix64F Pp) {
        InversionResult c;
        if (anyDiagonalInfinities(P)) {
            // Inflate variance
            final DenseMatrix64F Vp = matrix0;
            final DenseMatrix64F Vi = wrap(partials, bo + dimTrait + dimTrait * dimTrait, dimTrait, dimTrait);

            CommonOps.add(Vi, v, Vd, Vp);
            c = safeInvert(Vp, Pp, true);
        } else {
            // Deflate precision
            final DenseMatrix64F PPlusPd = matrix0;
            CommonOps.add(P, 1.0 / v, Pd, PPlusPd);

            final DenseMatrix64F PPlusPdInv = new DenseMatrix64F(dimTrait, dimTrait);
            safeInvert(PPlusPd, PPlusPdInv, false);

            CommonOps.mult(PPlusPdInv, P, Pp);
            CommonOps.mult(P, Pp, PPlusPdInv);
            CommonOps.add(P, -1, PPlusPdInv, Pp);
            c = safeDeterminant(Pp, false);
        }

        return c;
    }
}

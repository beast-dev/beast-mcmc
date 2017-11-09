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

//        assert precisionType == PrecisionType.FULL;
//
//        allocateStorage();
//
//        if (TIMING) {
//            times = new HashMap<String, Long>();
//        } else {
//            times = null;
//        }

        System.err.println("Trying SafeMultivariateIntegrator");
    }

//    @Override
//    public String getReport() {
//
//        StringBuilder sb = new StringBuilder();
//
//        if (TIMING) {
//            sb.append("\nTIMING:");
//            for (String key : times.keySet()) {
//                String value = String.format("%4.3e", (double) times.get(key));
//                sb.append("\n" + key + "\t\t" + value);
//            }
//            sb.append("\n");
//        }
//        return sb.toString();
//    }

    private static final boolean TIMING = false;

//    private final Map<String, Long> times;
//
//    private DenseMatrix64F matrix0;
//    private DenseMatrix64F matrix1;
//    private DenseMatrix64F matrix2;
//    private DenseMatrix64F matrix3;
//    private DenseMatrix64F matrix4;
//    private DenseMatrix64F matrix5;
//    private DenseMatrix64F matrix6;
//
//    private double[] vector0;

//    private void allocateStorage() {
//        inverseDiffusions = new double[dimTrait * dimTrait * diffusionCount];
//
//        vector0 = new double[dimTrait];
//        matrix0 = new DenseMatrix64F(dimTrait, dimTrait);
//        matrix1 = new DenseMatrix64F(dimTrait, dimTrait);
//        matrix2 = new DenseMatrix64F(dimTrait, dimTrait);
//        matrix3 = new DenseMatrix64F(dimTrait, dimTrait);
//        matrix4 = new DenseMatrix64F(dimTrait, dimTrait);
//        matrix5 = new DenseMatrix64F(dimTrait, dimTrait);
//        matrix6 = new DenseMatrix64F(dimTrait, dimTrait);
//    }

//    @Override
//    public void setDiffusionPrecision(int precisionIndex, final double[] matrix, double logDeterminant) {
//        super.setDiffusionPrecision(precisionIndex, matrix, logDeterminant);
//
//        assert (inverseDiffusions != null);
//
//        final int offset = dimTrait * dimTrait * precisionIndex;
//        DenseMatrix64F precision = wrap(diffusions, offset, dimTrait, dimTrait);
//        DenseMatrix64F variance = new DenseMatrix64F(dimTrait, dimTrait);
//        CommonOps.invert(precision, variance);
//        unwrap(variance, inverseDiffusions, offset);
//
//        if (DEBUG) {
//            System.err.println("At precision index: " + precisionIndex);
//            System.err.println("precision: " + precision);
//            System.err.println("variance : " + variance);
//        }
//    }

//    @Override
//    public boolean requireDataAugmentationForOuterProducts() {
//        return true;
//    }

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
//        // Read variance increments along descendent branches of k
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

        // Read variance increments along descendent branches of k
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
            final double lpip = Double.isInfinite(lpi) ?
                    1.0 / vi : lpi / (1.0 + lpi * vi);
            final double lpjp = Double.isInfinite(lpj) ?
                    1.0 / vj : lpj / (1.0 + lpj * vj);

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

            if (useVariancej) {

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

//            if (useVariance) {
//
//                final DenseMatrix64F Vip = matrix0;
//                final DenseMatrix64F Vjp = matrix1;
//
//                final DenseMatrix64F Vi = wrap(partials, ibo + dimTrait + dimTrait * dimTrait, dimTrait, dimTrait);
//                final DenseMatrix64F Vj = wrap(partials, jbo + dimTrait + dimTrait * dimTrait, dimTrait, dimTrait);
//
//                CommonOps.add(Vi, vi, Vd, Vip);
//                CommonOps.add(Vj, vj, Vd, Vjp);
//
//                ci = safeInvert(Vip, Pip, true);
//                cj = safeInvert(Vjp, Pjp, true);
//            } else {
//
//                final DenseMatrix64F PiPlusPd = matrix0;
//                final DenseMatrix64F PjPlusPd = matrix1;
//
//                CommonOps.add(Pi, 1.0 / vi, Pd, PiPlusPd);
//                CommonOps.add(Pj, 1.0 / vj, Pd, PjPlusPd);
//
//                final DenseMatrix64F PiPlusPdInv = new DenseMatrix64F(dimTrait, dimTrait);
//                final DenseMatrix64F PjPlusPdInv = new DenseMatrix64F(dimTrait, dimTrait);
//
//                safeInvert(PiPlusPd, PiPlusPdInv, false);
//                safeInvert(PjPlusPd, PjPlusPdInv, false);
//
//                CommonOps.mult(PiPlusPdInv, Pi, Pip);
//                CommonOps.mult(PjPlusPdInv, Pj, Pjp);
//
//                CommonOps.mult(Pi, Pip, PiPlusPdInv);
//                CommonOps.mult(Pj, Pjp, PjPlusPdInv);
//
//                CommonOps.add(Pi, -1, PiPlusPdInv, Pip);
//                CommonOps.add(Pj, -1, PjPlusPdInv, Pjp);
//
//                ci = safeDeterminant(Pip, false);
//                cj = safeDeterminant(Pjp, false);
//            }

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

//                final DenseMatrix64F Pk = new DenseMatrix64F(dimTrait, dimTrait);
            final DenseMatrix64F Pk = matrix4;
            CommonOps.add(Pip, Pjp, Pk);

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

//                final double[] tmp = new double[dimTrait];
            final double[] tmp = vector0;
            for (int g = 0; g < dimTrait; ++g) {
                double sum = 0.0;
                for (int h = 0; h < dimTrait; ++h) {
                    sum += Pip.unsafe_get(g, h) * partials[ibo + h];
                    sum += Pjp.unsafe_get(g, h) * partials[jbo + h];
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

            // TODO Use safeSolve just above

            if (TIMING) {
                endTime("peel4");
                startTime("peel5");
            }

            // C. Store precision
            partials[kbo + dimTrait + 2 * dimTrait * dimTrait] = lpk;

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
                System.err.print("\t\tmean j:");
                for (int e = 0; e < dimTrait; ++e) {
                    System.err.print(" " + partials[jbo + e]);
                }
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
                    final double ig = partials[ibo + g];
                    final double jg = partials[jbo + g];
                    final double kg = partials[kbo + g];

                    for (int h = 0; h < dimTrait; ++h) {
                        final double ih = partials[ibo + h];
                        final double jh = partials[jbo + h];
                        final double kh = partials[kbo + h];

                        SSi += ig * Pip.unsafe_get(g, h) * ih;
                        SSj += jg * Pjp.unsafe_get(g, h) * jh;
                        SSk += kg * Pk .unsafe_get(g, h) * kh;
                    }
                }

//                    final DenseMatrix64F Vt = new DenseMatrix64F(dimTrait, dimTrait);
//                final DenseMatrix64F Vt = matrix6;
//                CommonOps.add(Vip, Vjp, Vt);

//                if (DEBUG) {
//                    System.err.println("Vt: " + Vt);
//                }

                int dimensionChange = ci.getEffectiveDimension() + cj.getEffectiveDimension()
                        - ck.getEffectiveDimension();

//                    System.err.println(ci.getDeterminant());
//                    System.err.println(CommonOps.det(Vip));
//
//                    System.err.println(cj.getDeterminant());
//                    System.err.println(CommonOps.det(Vjp));
//
//                    System.err.println(1.0 / ck.getDeterminant());
//                    System.err.println(CommonOps.det(Vk));

                remainder += -dimensionChange * LOG_SQRT_2_PI - 0.5 *
//                            (Math.log(CommonOps.det(Vip)) + Math.log(CommonOps.det(Vjp)) - Math.log(CommonOps.det(Vk)))
                        (Math.log(ci.getDeterminant()) + Math.log(cj.getDeterminant()) + Math.log(ck.getDeterminant()))
                        - 0.5 * (SSi + SSj - SSk);

                // TODO Can get SSi + SSj - SSk from inner product w.r.t Pt (see outer-products below)?

                if (DEBUG) {
                    System.err.println("\t\t\tSSi = " + (SSi));
                    System.err.println("\t\t\tSSj = " + (SSj));
                    System.err.println("\t\t\tSSk = " + (SSk));
                    System.err.println("\t\t\tdeti = " + Math.log(ci.getDeterminant()));
                    System.err.println("\t\t\tdetj = " + Math.log(ci.getDeterminant()));
                    System.err.println("\t\t\tdetk = " + Math.log(ci.getDeterminant()));
                    System.err.println("\t\tremainder: " + remainder);
//                        System.exit(-1);
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

            final DenseMatrix64F Proot = wrap(partials, rootOffset + dimTrait, dimTrait, dimTrait);
            final DenseMatrix64F Pprior = wrap(partials, priorOffset + dimTrait, dimTrait, dimTrait);

//            final DenseMatrix64F Vroot = wrap(partials, rootOffset + dimTrait + dimTrait * dimTrait, dimTrait, dimTrait);
//            final DenseMatrix64F Vprior = wrap(partials, priorOffset + dimTrait + dimTrait * dimTrait, dimTrait, dimTrait);

            // TODO Block below is for the conjugate prior ONLY
            {
//                final DenseMatrix64F Vtmp = new DenseMatrix64F(dimTrait, dimTrait);
//                CommonOps.mult(Vd, Vprior, Vtmp);
//                Vprior.set(Vtmp);

                final DenseMatrix64F Ptmp = new DenseMatrix64F(dimTrait, dimTrait);
                CommonOps.mult(Pd, Pprior, Ptmp);
                Pprior.set(Ptmp); // TODO What does this do?
            }

            final DenseMatrix64F Vtotal = new DenseMatrix64F(dimTrait, dimTrait);
//            CommonOps.add(Vroot, Vprior, Vtotal);

            final DenseMatrix64F Ptotal = new DenseMatrix64F(dimTrait, dimTrait);
            CommonOps.invert(Vtotal, Ptotal);  // TODO Can return determinant at same time to avoid extra QR decomp

            final DenseMatrix64F tmp1 = new DenseMatrix64F(dimTrait, dimTrait);
            final DenseMatrix64F tmp2 = new DenseMatrix64F(dimTrait, dimTrait);
            CommonOps.add(Proot, Pprior, Ptotal);
            CommonOps.invert(Ptotal, Vtotal);
            CommonOps.mult(Vtotal, Proot, tmp1);
            CommonOps.mult(Proot, tmp1, tmp2);
            CommonOps.add(Proot, -1.0, tmp2, Ptotal);

            double SS = 0;
            for (int g = 0; g < dimTrait; ++g) {
                final double gDifference = partials[rootOffset + g] - partials[priorOffset + g];

                for (int h = 0; h < dimTrait; ++h) {
                    final double hDifference = partials[rootOffset + h] - partials[priorOffset + h];

                    SS += gDifference * Ptotal.unsafe_get(g, h) * hDifference;
                }
            }

            final double logLike = -dimTrait * LOG_SQRT_2_PI
//                    - 0.5 * Math.log(CommonOps.det(Vtotal))
                    + 0.5 * Math.log(CommonOps.det(Ptotal))
                    - 0.5 * SS;

            final double remainder = remainders[rootBufferIndex * numTraits + trait];
            logLikelihoods[trait] = logLike + remainder;

            if (incrementOuterProducts) {

                assert false : "Should not get here";

//                int opo = dimTrait * dimTrait * trait;
//                int opd = precisionOffset;
//
//                double rootScalar = partials[rootOffset + dimTrait + 2 * dimTrait * dimTrait];
//                final double priorScalar = partials[priorOffset + dimTrait];
//
//                if (!Double.isInfinite(priorScalar)) {
//                    rootScalar = rootScalar * priorScalar / (rootScalar + priorScalar);
//                }
//
//                for (int g = 0; g < dimTrait; ++g) {
//                    final double gDifference = partials[rootOffset + g] - partials[priorOffset + g];
//
//                    for (int h = 0; h < dimTrait; ++h) {
//                        final double hDifference = partials[rootOffset + h] - partials[priorOffset + h];
//
//                        outerProducts[opo] += gDifference * hDifference
////                                    * Ptotal.unsafe_get(g, h) / diffusions[opd];
//                                * rootScalar;
//                        ++opo;
//                        ++opd;
//                    }
//                }
//
//                degreesOfFreedom[trait] += 1; // incremenent degrees-of-freedom
            }

            if (DEBUG) {
                System.err.print("mean:");
                for (int g = 0; g < dimTrait; ++g) {
                    System.err.print(" " + partials[rootOffset + g]);
                }
                System.err.println("");
                System.err.println("Proot: " + Proot);
//                System.err.println("Vroot: " + Vroot);
                System.err.println("Pprior: " + Pprior);
//                System.err.println("Vprior: " + Vprior);
                System.err.println("Ptotal: " + Ptotal);
//                System.err.println("Vtotal: " + Vtotal);
//                    System.err.println("prec: " + partials[rootOffset + dimTrait]);
                System.err.println("\t" + logLike + " " + (logLike + remainder));

                if (incrementOuterProducts) {
                    System.err.println("Outer-products:" + wrap(outerProducts, dimTrait * dimTrait * trait, dimTrait, dimTrait));
                }
            }

            rootOffset += dimPartialForTrait;
            priorOffset += dimPartialForTrait;
        }

        if (DEBUG) {
            System.err.println("End");
        }
    }
}

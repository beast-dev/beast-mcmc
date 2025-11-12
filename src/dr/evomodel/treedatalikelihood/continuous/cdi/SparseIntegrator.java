/*
 * SparseIntegrator.java
 *
 * Copyright © 2002-2025 the BEAST Development Team
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

import java.util.HashMap;
import java.util.Map;

import static dr.math.matrixAlgebra.missingData.MissingOps.wrap;

/**
 * @author Marc A. Suchard
 * @author Philippe Lemey
 */

public class SparseIntegrator extends ContinuousDiffusionIntegrator.Basic {

    private static final boolean DEBUG = true;
    private static final boolean TIMING = true;

    private double[] sumOfSquares;

    private final Map<String, Long> times;

    public SparseIntegrator(PrecisionType precisionType, int numTraits, int dimTrait, int dimProcess,
                            int bufferCount, int diffusionCount) {
        super(precisionType, numTraits, dimTrait, dimProcess, bufferCount, diffusionCount);

        assert precisionType == PrecisionType.FULL;

        allocateStorage();

        if (TIMING) {
            times = new HashMap<>();
        } else {
            times = null;
        }
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
        // Determine buffer offsets
        int kbo = dimPartial * kBuffer;
        int ibo = dimPartial * iBuffer;
        int jbo = dimPartial * jBuffer;

        // Determine matrix offsets
        final int imo = iMatrix; //TODO: just use iMatrix * jMatrix? (also, why do we need these?)
        final int jmo = jMatrix;

        // Read variance increments along descendant branches of k
        final double vi = branchLengths[imo];
        final double vj = branchLengths[jmo];

        if (DEBUG) {
            System.err.println("i:");
            System.err.println("\tvar : " + branchLengths[imo]);
        }

        // For each trait // TODO in parallel
        for (int trait = 0; trait < numTraits; ++trait) {

            // Increase variance along the branches i -> k and j -> k

            // A. Get current precision of i and j
            final double pi = partials[ibo + dimTrait];
            final double pj = partials[jbo + dimTrait];

            // B. Integrate along branch using two matrix inversions
            final double pip = Double.isInfinite(pi) ?
                    1.0 / vi : pi / (1.0 + pi * vi);
            final double pjp = Double.isInfinite(pj) ?
                    1.0 / vj : pj / (1.0 + pj * vj);

            // Compute partial mean and precision at node k

            // A. Partial precision scalar
            final double pk = Double.isInfinite(pip + pjp) ? Double.POSITIVE_INFINITY : pip + pjp;
            final double pTmp = Double.isInfinite(pip) ? 1.0 : 1.0 + pjp / pip;
            final double pWeight = Double.isInfinite(pTmp) ? 0.0 : 1.0 / pTmp;

            // B. Partial mean
            if (INLINE) {
                // For each dimension // TODO in parallel
                for (int g = 0; g < dimTrait; ++g) {
//                        partials[kbo + g] = (pip * partials[ibo + g] + pjp * partials[jbo + g]) / pk;
                    partials[kbo + g] = pWeight * partials[ibo + g] + (1.0 - pWeight) * partials[jbo + g];
                }
            } else {
                updateMean(partials, kbo, ibo, jbo, pip, pjp, pk, dimTrait);
            }

            // C. Store precision
            partials[kbo + dimTrait] = pk;

            if (DEBUG) {
                System.err.println("\ttrait: " + trait);
                //System.err.println("\t\tPrec: " + pi);
                System.err.print("\t\tmean i:");
                for (int e = 0; e < dimTrait; ++e) {
                    System.err.print(" " + partials[ibo + e]);
                }
                System.err.println(" prec i: " + pi);
                System.err.print("\t\tmean j:");
                for (int e = 0; e < dimTrait; ++e) {
                    System.err.print(" " + partials[jbo + e]);
                }
                System.err.println(" prec j: " + pj);

                if (pj == 0.0) { System.exit(-1); }
                System.err.print("\t\tmean k:");
                for (int e = 0; e < dimTrait; ++e) {
                    System.err.print(" " + partials[kbo + e]);
                }
                System.err.println(" prec k: " + pk);
                System.err.println("");
            }

            // Computer remainder at node k
            double remainder = 0.0;
            if (computeRemainders && pi != 0 && pj != 0) {

                // TODO Suspect this is very inefficient, since SSi and SSj were already computed when k = i or j

                remainder += computeRemainder(trait, ibo, jbo, kbo, pip, pjp, pk, incrementOuterProducts);

//                final double remainderPrecision = pip * pjp / pk;
//
//                // Inner products
//                double SSk = 0;
//                double SSj = 0;
//                double SSi = 0;
//
//                int pob = precisionOffset;
//
//                // vector-matrix-vector TODO in parallel
//                for (int g = 0; g < dimTrait; ++g) {
//                    final double ig = partials[ibo + g];
//                    final double jg = partials[jbo + g];
//                    final double kg = partials[kbo + g];
//
//                    for (int h = 0; h < dimTrait; ++h) {
//                        final double ih = partials[ibo + h];
//                        final double jh = partials[jbo + h];
//                        final double kh = partials[kbo + h];
//
//                        final double diffusion = diffusions[pob]; // element [g][h]
//
//                        SSi += ig * diffusion * ih;
//                        SSj += jg * diffusion * jh;
//                        SSk += kg * diffusion * kh;
//
//                        ++pob;
//                    }
//                }
//
//                remainder += -dimTrait * LOG_SQRT_2_PI // TODO Can move some calculation outside the loop
//                        + 0.5 * (dimTrait * Math.log(remainderPrecision) + precisionLogDet)
//                        - 0.5 * (pip * SSi + pjp * SSj - pk * SSk);
//
//                if (DEBUG) {
//                    System.err.println("\t\t\tpk: " + pk);
//                    System.err.println("\t\t\tSSi = " + (pip * SSi));
//                    System.err.println("\t\t\tSSj = " + (pjp * SSj));
//                    System.err.println("\t\t\tSSk = " + (pk * SSk));
//                }
//
//
//                if (DEBUG) {
//                    System.err.println("\t\tremainder: " + remainder);
//                }
//
//                if (incrementOuterProducts) {
//                    int opo = dimTrait * dimTrait * trait;
//
////                        final double remainderPrecision = pip * pjp / (pip + pjp);
//
//                    if (DEBUG) {
//                        System.err.println("pip: " + pip);
//                        System.err.println("pjp: " + pjp);
//                        System.err.println("sum: " + (pip + pjp));
//                        System.err.println("op prec: " + remainderPrecision);
//                    }
//
//                    for (int g = 0; g < dimTrait; ++g) {
//                        final double ig = partials[ibo + g];
//                        final double jg = partials[jbo + g];
//
//                        for (int h = 0; h < dimTrait; ++h) {
//                            final double ih = partials[ibo + h];
//                            final double jh = partials[jbo + h];
//
//                            outerProducts[opo] += (ig - jg) * (ih - jh) * remainderPrecision;
//                            ++opo;
//                        }
//                    }
//
//                    if (DEBUG) {
//                        System.err.println("Outer-products:" + wrap(outerProducts, dimTrait * dimTrait * trait, dimTrait, dimTrait));
//                    }
//
//                    degreesOfFreedom[trait] += 1; // increment degrees-of-freedom
//                }
            } // End if remainder

            // Accumulate remainder up tree and store

            remainders[kBuffer * numTraits + trait] = remainder
                    + remainders[iBuffer * numTraits + trait] + remainders[jBuffer * numTraits + trait];

            // Get ready for next trait
            kbo += dimPartialForTrait;
            ibo += dimPartialForTrait;
            jbo += dimPartialForTrait;

        }
    }

    private void allocateStorage() {
        sumOfSquares = new double[numTraits * bufferCount];
    }

    private double computeRemainder(
            final int trait,
            final int ibo,
                             final int jbo,
                             final int kbo,
                             final double pip,
                             final double pjp,
                             final double pk,
            final boolean incrementOuterProducts) {

        // TODO Suspect this is very inefficient, since SSi and SSj were already computed when k = i or j

        final double remainderPrecision = pip * pjp / pk;

        double remainder = 0;

        // Inner products
        double SSk = 0;
        double SSj = 0;
        double SSi = 0;

        int pob = precisionOffset;

        // vector-matrix-vector TODO in parallel
        for (int g = 0; g < dimTrait; ++g) {
            final double ig = partials[ibo + g];
            final double jg = partials[jbo + g];
            final double kg = partials[kbo + g];

            for (int h = 0; h < dimTrait; ++h) {
                final double ih = partials[ibo + h];
                final double jh = partials[jbo + h];
                final double kh = partials[kbo + h];

                final double diffusion = diffusions[pob]; // element [g][h]

                SSi += ig * diffusion * ih;
                SSj += jg * diffusion * jh;
                SSk += kg * diffusion * kh;

                ++pob;
            }
        }

        // TODO Test post-order compute of SSs
        sumOfSquares[kbo + trait] = SSk;

        remainder += -dimTrait * LOG_SQRT_2_PI // TODO Can move some calculation outside the loop
                + 0.5 * (dimTrait * Math.log(remainderPrecision) + precisionLogDet)
                - 0.5 * (pip * SSi + pjp * SSj - pk * SSk);

        if (DEBUG) {
            System.err.println("\t\t\tpk: " + pk);
            System.err.println("\t\t\tSSi = " + (pip * SSi));
            System.err.println("\t\t\tSSj = " + (pjp * SSj));
            System.err.println("\t\t\tSSk = " + (pk * SSk));
        }


        if (DEBUG) {
            System.err.println("\t\tremainder: " + remainder);
        }

        if (incrementOuterProducts) {
            int opo = dimTrait * dimTrait * trait;

//                        final double remainderPrecision = pip * pjp / (pip + pjp);

            if (DEBUG) {
                System.err.println("pip: " + pip);
                System.err.println("pjp: " + pjp);
                System.err.println("sum: " + (pip + pjp));
                System.err.println("op prec: " + remainderPrecision);
            }

            for (int g = 0; g < dimTrait; ++g) {
                final double ig = partials[ibo + g];
                final double jg = partials[jbo + g];

                for (int h = 0; h < dimTrait; ++h) {
                    final double ih = partials[ibo + h];
                    final double jh = partials[jbo + h];

                    outerProducts[opo] += (ig - jg) * (ih - jh) * remainderPrecision;
                    ++opo;
                }
            }

            if (DEBUG) {
                System.err.println("Outer-products:" + wrap(outerProducts, dimTrait * dimTrait * trait, dimTrait, dimTrait));
            }

            degreesOfFreedom[trait] += 1; // increment degrees-of-freedom
        }

        return remainder;
    }
}

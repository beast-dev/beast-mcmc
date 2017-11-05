/*
 * ContinuousDiffusionIntegrator.java
 *
 * Copyright (c) 2002-2016 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
 */

package dr.evomodel.treedatalikelihood.continuous.cdi;

import dr.math.matrixAlgebra.WrappedVector;
import dr.xml.Reportable;

import static dr.math.matrixAlgebra.missingData.MissingOps.*;

/**
 * @author Marc A. Suchard
 */
public interface ContinuousDiffusionIntegrator extends Reportable {
    int OPERATION_TUPLE_SIZE = 5;
    int NONE = -1;

    double LOG_SQRT_2_PI = 0.5 * Math.log(2 * Math.PI);

    void finalize() throws Throwable;

    void setPostOrderPartial(int bufferIndex, final double[] partial);

    void getPostOrderPartial(int bufferIndex, final double[] partial);

    double getBranchMatrices(int bufferIndex, final double[] precision, final double[] displacement);

    void setPreOrderPartial(int bufferIndex, final double[] partial);

    void getPreOrderPartial(int bufferIndex, final double[] partial);

    void setWishartStatistics(final int[] degreesOfFreedom, final double[] outerProducts);

    void getWishartStatistics(final int[] degreesOfFreedom, final double[] outerProducts);

    void setDiffusionPrecision(int diffusionIndex, final double[] matrix, double logDeterminant);

    void updatePostOrderPartials(final int[] operations, int operationCount, boolean incrementOuterProducts);

    void updatePreOrderPartials(final int[] operations, int operationCount);

    InstanceDetails getDetails();

    void updateBrownianDiffusionMatrices(int precisionIndex, final int[] probabilityIndices,
                                         final double[] edgeLengths, final double[] driftRates,
                                         int updateCount);

//    void updateOrnsteinUhlenbeckMatrices(int precisionIndex, final int[] probabilityIndices,
//                                         final double[] edgeLengths,
//                                         int updateCount);

    void calculateRootLogLikelihood(int rootBufferIndex, int priorBufferIndex, double[] logLike,
                                    boolean incrementOuterProducts);

//    int getPartialBufferCount();

//    int getMatrixBufferCount();

//    boolean requireDataAugmentationForOuterProducts();

    // TODO Only send a list of operations
    void updatePreOrderPartial(int kp, int ip, int im, int jp, int jm);

    class Basic implements ContinuousDiffusionIntegrator {

//        private int instance = -1;
        private InstanceDetails details = new InstanceDetails();

//        private final PrecisionType precisionType;
        final int numTraits;
        final int dimTrait;
        final int bufferCount;
        final int diffusionCount;

        final int dimMatrix;
        final int dimPartialForTrait;
        final int dimPartial;

        @Override
        public String getReport() {
            return "";
        }

        public Basic(
                final PrecisionType precisionType,
                final int numTraits,
                final int dimTrait,
                final int bufferCount,
                final int diffusionCount
        ) {
            assert(numTraits > 0);
            assert(dimTrait > 0);
            assert(bufferCount > 0);
            assert(diffusionCount > 0);

//            assert (precisionType == PrecisionType.SCALAR);

//            this.precisionType = precisionType;
            this.numTraits = numTraits;
            this.dimTrait = dimTrait;
            this.bufferCount = bufferCount;
            this.diffusionCount = diffusionCount;

            this.dimMatrix = precisionType.getMatrixLength(dimTrait);
            this.dimPartialForTrait = dimTrait + dimMatrix;
            this.dimPartial = numTraits * dimPartialForTrait;

            if (DEBUG) {
                System.err.println("numTraits: " + numTraits);
                System.err.println("dimTrait: " + dimTrait);
                System.err.println("dimMatrix: " + dimMatrix);
                System.err.println("dimPartialForTrait: " + dimPartialForTrait);
                System.err.println("dimPartial: " + dimPartial);
            }

            allocateStorage();
        }

//        @Override
//        public int getPartialBufferCount() { return  bufferCount; }
//
//        @Override
//        public boolean requireDataAugmentationForOuterProducts() {
//            return false;
//        }

        @Override
        public void finalize() throws Throwable {
            super.finalize();
        }

        @Override
        public void setPostOrderPartial(int bufferIndex, final double[] partial) {
            assert(partial.length == dimPartial);
            assert(partials != null);

            System.arraycopy(partial, 0, partials, dimPartial * bufferIndex, dimPartial);
        }

        @Override
        public void getPostOrderPartial(int bufferIndex, final double[] partial) {
            assert(partials != null);
            assert(partial.length >= dimPartial);

            System.arraycopy(partials, dimPartial * bufferIndex, partial, 0, dimPartial);
        }
        
        @Override
        public double getBranchMatrices(int bufferIndex, double[] precision, double[] displacement) {
            return 1.0 / branchLengths[bufferIndex * dimMatrix];
        }

        @Override
        public void setPreOrderPartial(int bufferIndex, final double[] partial) {
            assert(partial.length >= dimPartial);
            assert(prePartials != null);

            System.arraycopy(partial, 0, prePartials, dimPartial * bufferIndex, dimPartial);
        }

        @Override
        public void getPreOrderPartial(int bufferIndex, final double[] partial) {
            assert(partial.length >= dimPartial);
            assert(prePartials != null);

            System.arraycopy(prePartials, dimPartial * bufferIndex, partial, 0, dimPartial);
        }

        @Override
        public void setWishartStatistics(final int[] degreesOfFreedom, final double[] outerProducts) {
            assert(degreesOfFreedom.length == numTraits);
            assert(outerProducts.length == dimTrait * dimTrait * numTraits);

            System.arraycopy(degreesOfFreedom, 0, this.degreesOfFreedom, 0, numTraits);
            System.arraycopy(outerProducts, 0, this.outerProducts, 0, dimTrait * dimTrait * numTraits);
        }

        @Override
        public void getWishartStatistics(final int[] degreesOfFreedom, final double[] outerProducts) {
            assert(degreesOfFreedom.length == numTraits);
            assert(outerProducts.length == dimTrait * dimTrait * numTraits);

            System.arraycopy(this.degreesOfFreedom, 0, degreesOfFreedom, 0, numTraits);
            System.arraycopy(this.outerProducts, 0, outerProducts, 0, dimTrait * dimTrait * numTraits);
        }

        @Override
        public void setDiffusionPrecision(int precisionIndex, final double[] matrix, double logDeterminant) {
            assert(matrix.length == dimTrait * dimTrait);
            assert(diffusions != null);
            assert(determinants != null);

            System.arraycopy(matrix, 0, diffusions, dimTrait * dimTrait * precisionIndex, dimTrait * dimTrait);
            determinants[precisionIndex] = logDeterminant;
        }

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

            // TODO For each trait in parallel
            for (int trait = 0; trait < numTraits; ++trait) {

                double SS = 0;
                int pob = precisionOffset;

                double rootScalar = partials[rootOffset + dimTrait];
                final double priorScalar = partials[priorOffset + dimTrait];

                if (!Double.isInfinite(priorScalar)) {
                    rootScalar = rootScalar * priorScalar / (rootScalar + priorScalar);
                }

                for (int g = 0; g < dimTrait; ++g) {
                    final double gDifference = partials[rootOffset + g] - partials[priorOffset + g];

                    for (int h = 0; h < dimTrait; ++h) {
                        final double hDifference = partials[rootOffset + h] - partials[priorOffset + h];

                        SS += gDifference * diffusions[pob] * hDifference;
                        ++pob;
                    }
                }

                final double logLike = -dimTrait * LOG_SQRT_2_PI
                        + 0.5 * (dimTrait * Math.log(rootScalar) + precisionLogDet)
                        - 0.5 * rootScalar * SS;
                final double remainder = remainders[rootBufferIndex * numTraits + trait];

                logLikelihoods[trait] = logLike + remainder;

                if (incrementOuterProducts) {
                    int opo = dimTrait * dimTrait * trait;

                    for (int g = 0; g < dimTrait; ++g) {
                        final double gDifference = partials[rootOffset + g] - partials[priorOffset + g];

                        for (int h = 0; h < dimTrait; ++h) {
                            final double hDifference = partials[rootOffset + h] - partials[priorOffset + h];

                            outerProducts[opo] += gDifference * hDifference * rootScalar;
                            ++opo;
                        }
                    }

                    degreesOfFreedom[trait] += 1; // incremenent degrees-of-freedom
                }

                if (DEBUG) {
                    System.err.print("mean:");
                    for (int g = 0; g < dimTrait; ++g) {
                        System.err.print(" " + partials[rootOffset + g]);
                    }
                    System.err.println("");
                    System.err.println("prec: " + partials[rootOffset + dimTrait]);
                    System.err.println("rootScalar: " + rootScalar);
                    System.err.println("\t" + logLike + " " + (logLike + remainder));
                    if (incrementOuterProducts) {
                        System.err.println("Outer-products:" + wrap(outerProducts, dimTrait * dimTrait * trait, dimTrait, dimTrait));
                    }
                    System.err.println("");
                }

                rootOffset += dimPartialForTrait;
                priorOffset += dimPartialForTrait;
            }

            if (DEBUG) {
                System.err.println("End");
            }
        }

        @Override
        public void updatePreOrderPartials(final int[] operations, int operationCount) {

            if (DEBUG) {
                System.err.println("Pre-order operations:");
            }


            throw new RuntimeException("Not yet implemented");
        }

        @Override
        public void updatePostOrderPartials(final int[] operations, int operationCount, boolean incrementOuterProducts) {

            if (DEBUG) {
                System.err.println("Post-order operations:");
            }

            int offset = 0;
            for (int op = 0; op < operationCount; ++op) {

                if (DEBUG) {
                    System.err.println("\t" + getOperationString(operations, offset));
                }

                updatePartial(
                        operations[offset    ],
                        operations[offset + 1],
                        operations[offset + 2],
                        operations[offset + 3],
                        operations[offset + 4],
                        incrementOuterProducts
                );

                offset += ContinuousDiffusionIntegrator.OPERATION_TUPLE_SIZE;
            }

            if (DEBUG) {
                System.err.println("End");
                System.err.println("");
            }
        }

        @Override
        public void updateBrownianDiffusionMatrices(int precisionIndex, final int[] probabilityIndices,
                                                    final double[] edgeLengths, final double[] driftRates,
                                                    int updateCount) {

            if (DEBUG) {
                System.err.println("Matrices (basic):");
            }

            for (int up = 0; up < updateCount; ++up) {

                if (DEBUG) {
                    System.err.println("\t" + probabilityIndices[up] + " <- " + edgeLengths[up]);
                }

                // TODO Currently only written for SCALAR model
                branchLengths[dimMatrix * probabilityIndices[up]] = edgeLengths[up];  // TODO Remove dimMatrix
            }

            precisionOffset = dimTrait * dimTrait * precisionIndex;
            precisionLogDet = determinants[precisionIndex];
        }

        @Override
        public InstanceDetails getDetails() {
            return details;
        }

        // Internal storage
        double[] partials;

        double[] branchLengths;
        double[] variances;
        double[] precisions;

        double[] remainders;
        double[] diffusions;
        double[] determinants;
        int[] degreesOfFreedom;
        double[] outerProducts;
        double[] prePartials;

        // Set during updateDiffusionMatrices() and used in updatePartials()
        int precisionOffset;
        double precisionLogDet;

        static final boolean INLINE = true;

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
            final int imo = dimMatrix * iMatrix;
            final int jmo = dimMatrix * jMatrix;

            // Read variance increments along descendent branches of k
            final double vi = branchLengths[imo];
            final double vj = branchLengths[jmo];
//
//            final DenseMatrix64F Vd = wrap(inverseDiffusions, precisionOffset, dimTrait, dimTrait);
//
            if (DEBUG) {
                System.err.println("updatePreOrderPartial for node " + iBuffer);
                System.err.println("\tvi: " + vi + " vj: " + vj);
            }

            // For each trait // TODO in parallel
            for (int trait = 0; trait < numTraits; ++trait) {

                // A. Get current precision of k and j
                final double pk = prePartials[kbo + dimTrait];
                final double pj = partials[jbo + dimTrait];

//                final DenseMatrix64F Pk = wrap(prePartials, kbo + dimTrait, dimTrait, dimTrait);
//    //                final DenseMatrix64F Pj = wrap(partials, jbo + dimTrait, dimTrait, dimTrait);
//
//    //                final DenseMatrix64F Vk = wrap(prePartials, kbo + dimTrait + dimTrait * dimTrait, dimTrait, dimTrait);
//                final DenseMatrix64F Vj = wrap(partials, jbo + dimTrait + dimTrait * dimTrait, dimTrait, dimTrait);
//

                // B. Inflate variance along sibling branch using matrix inversion
                final double pjp = Double.isInfinite(pj) ?
                        1.0 / vj : pj / (1.0 + pj * vj);

//    //                final DenseMatrix64F Vjp = new DenseMatrix64F(dimTrait, dimTrait);
//                final DenseMatrix64F Vjp = matrix0;
//                CommonOps.add(Vj, vj, Vd, Vjp);
//
//    //                final DenseMatrix64F Pjp = new DenseMatrix64F(dimTrait, dimTrait);
//                final DenseMatrix64F Pjp = matrix1;
//                InversionResult cj = safeInvert(Vjp, Pjp, false);
//
//    //                final DenseMatrix64F Pip = new DenseMatrix64F(dimTrait, dimTrait);
//                final DenseMatrix64F Pip = matrix2;
//                CommonOps.add(Pk, Pjp, Pip);

                final double pip = pjp + pk;
//
//    //                final DenseMatrix64F Vip = new DenseMatrix64F(dimTrait, dimTrait);
//                final DenseMatrix64F Vip = matrix3;
//                InversionResult cip = safeInvert(Pip, Vip, false);
//
//                // C. Compute prePartial mean
//    //                final double[] tmp = new double[dimTrait];
//                final double[] tmp = vector0;
//                for (int g = 0; g < dimTrait; ++g) {
//                    double sum = 0.0;
//                    for (int h = 0; h < dimTrait; ++h) {
//                        sum += Pk.unsafe_get(g, h) * prePartials[kbo + h]; // Read parent
//                        sum += Pjp.unsafe_get(g, h) * partials[jbo + h];   // Read sibling
//                    }
//                    tmp[g] = sum;
//                }
//                for (int g = 0; g < dimTrait; ++g) {
//                    double sum = 0.0;
//                    for (int h = 0; h < dimTrait; ++h) {
//                        sum += Vip.unsafe_get(g, h) * tmp[h];
//                    }
//                    prePartials[ibo + g] = sum; // Write node
//                }

                for (int g = 0; g < dimTrait; ++g) {
                    prePartials[ibo + g] = (pk * prePartials[kbo + g] + pjp * partials[jbo + g]) / pip;
                }

                // C. Inflate variance along node branch

                final double pi  = Double.isInfinite(pip) ?
                        1.0 / vi : pip / (1.0 + pip * vi);

//                final DenseMatrix64F Vi = Vip;
//                CommonOps.add(vi, Vd, Vip, Vi);
//
//    //                final DenseMatrix64F Pi = new DenseMatrix64F(dimTrait, dimTrait);
//                final DenseMatrix64F Pi = matrix4;
//                InversionResult ci = safeInvert(Vi, Pi, false);
//
//                // X. Store precision results for node

                prePartials[ibo + dimTrait] = pi;

//                unwrap(Pi, prePartials, ibo + dimTrait);
//                unwrap(Vi, prePartials, ibo + dimTrait + dimTrait * dimTrait);
//
                if (DEBUG) {
                    System.err.println("trait: " + trait);
                    System.err.println("pM: " + new WrappedVector.Raw(prePartials, kbo, dimTrait));
                    System.err.println("pk: " + pk);
                    System.err.println("sM: " + new WrappedVector.Raw(partials, jbo, dimTrait));
                    System.err.println("sV: " + vj);
//                    System.err.println("sVp: " + Vjp);
                    System.err.println("sPp: " + pjp);
                    System.err.println("Pip: " + pip);
                    System.err.println("cM: " + new WrappedVector.Raw(prePartials, ibo, dimTrait));
                    System.err.println("cV: " + pi);
                }
//
                // Get ready for next trait
                kbo += dimPartialForTrait;
                ibo += dimPartialForTrait;
                jbo += dimPartialForTrait;
            }
//            throw new RuntimeException("Not yet implemented");
        }

        protected void updatePartial(
                final int kBuffer,
                final int iBuffer,
                final int iMatrix,
                final int jBuffer,
                final int jMatrix,
                final boolean incrementOuterProducts
        ) {
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
                final double pk = pip + pjp;

                // B. Partial mean
                if (INLINE) {
                    // For each dimension // TODO in parallel
                    for (int g = 0; g < dimTrait; ++g) {
                        partials[kbo + g] = (pip * partials[ibo + g] + pjp * partials[jbo + g]) / pk;
                    }
                } else {
                    updateMean(partials, kbo, ibo, jbo, pip, pjp, pk, dimTrait);
                }

                // C. Store precision
                partials[kbo + dimTrait] = pk;

                if (DEBUG) {
                    System.err.println("\ttrait: " + trait);
                    //System.err.println("\t\tprec: " + pi);
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
                if (pi != 0 && pj != 0) {

                    // TODO Suspect this is very inefficient, since SSi and SSj were already computed when k = i or j

                    final double remainderPrecision = pip * pjp / pk;

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

                        degreesOfFreedom[trait] += 1; // incremenent degrees-of-freedom
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
        }

        private static void updateMean(final double[] partials,
                                       final int kob,
                                       final int iob,
                                       final int job,
                                       final double pip,
                                       final double pjp,
                                       final double pk,
                                       final int dimTrait) {
            for (int g = 0; g < dimTrait; ++g) {
                partials[kob + g] = (pip * partials[iob + g] + pjp * partials[job + g]) / pk;
            }
        }

        private void allocateStorage() {
            partials = new double[dimPartial * bufferCount];
            branchLengths = new double[dimMatrix * bufferCount]; // TODO Should be just bufferCount
//            variances = new double[dimMatrix * bufferCount]; // TODO Should be dimTrait * dimTrait
            remainders = new double[numTraits * bufferCount];

            diffusions = new double[dimTrait * dimTrait * diffusionCount];
            determinants = new double[diffusionCount];

            degreesOfFreedom = new int[numTraits];
            outerProducts = new double[dimTrait * dimTrait * numTraits];

            prePartials = new double[dimPartial * bufferCount];
        }

        private String getOperationString(final int[] operations, final int offset) {
            StringBuilder sb = new StringBuilder("op:");
            for (int i = 0; i < ContinuousDiffusionIntegrator.OPERATION_TUPLE_SIZE; ++i) {
                sb.append(" ").append(operations[offset + i]);
            }
            return sb.toString();
        }

        private static boolean DEBUG = false;
    }

}

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

/**
 * @author Marc A. Suchard
 */
public interface ContinuousDiffusionIntegrator {
    int OPERATION_TUPLE_SIZE = 5;
    int NONE = -1;

    double LOG_SQRT_2_PI = 0.5 * Math.log(2 * Math.PI);

    void finalize() throws Throwable;

//    void setPartialMean(int bufferIndex, final double[] mean);

//    void getPartialMean(int bufferIndex, final double[] mean);

//    void setPartialPrecision(int bufferIndex, final double[] precision);

//    void getPartialPrecision(int bufferIndex, final double[] precision);

    void setPartial(int bufferIndex, final double[] partial);

    void getPartial(int bufferIndex, final double[] partial);

    void setWishartStatistics(final int[] degreesOfFreedom, final double[] outerProducts);

    void getWishartStatistics(final int[] degreesOfFreedom, final double[] outerProducts);

    void setDiffusionPrecision(int diffusionIndex, final double[] matrix, double logDeterminant);

    void updatePartials(final int[] operations, int operationCount, boolean incrementOuterProducts);

    InstanceDetails getDetails();

    void updateDiffusionMatrices(int precisionIndex, final int[] probabilityIndices, final double[] edgeLengths,
                                 int updateCount);

    void calculateRootLogLikelihood(int rootBufferIndex, int priorBufferIndex, double[] logLike,
                                    boolean incrementOuterProducts);

//    abstract class AbstractBase implements ContinuousDiffusionIntegrator {
//
//        private InstanceDetails details = new InstanceDetails();
//
//        @Override
//        public void finalize() throws Throwable {
//            super.finalize();
//        }
//
//        @Override
//        public void setPartialMean(int bufferIndex, double[] mean) {
//
//        }
//
//        @Override
//        public void getPartialMean(int bufferIndex, double[] mean) {
//
//        }
//
//        @Override
//        public void setPartialPrecision(int bufferIndex, double[] precision) {
//
//        }
//
//        @Override
//        public void getPartialPrecision(int bufferIndex, double[] precision) {
//
//        }
//
//
//
//        @Override
//        public abstract void setPartial(int bufferIndex, double[] partial);
//
//        @Override
//        public abstract void getPartial(int bufferIndex, double[] partial);
//
//        @Override
//        public abstract void setDiffusionPrecision(int diffusionIndex, double[] matrix);
//
//        @Override
//        public abstract void updatePartials(int[] operations, int operationCount);
//
//    }

    class Basic implements ContinuousDiffusionIntegrator {

        private int instance = -1;
        private InstanceDetails details = new InstanceDetails();

        private final PrecisionType precisionType;
        private final int numTraits;
        private final int dimTrait;
        private final int bufferCount;
        private final int diffusionCount;

        private final int dimMatrix;
        private final int dimPartialForTrait;
        private final int dimPartial;

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

            this.precisionType = precisionType;
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

        @Override
        public void finalize() throws Throwable {
            super.finalize();
        }

        @Override
        public void setPartial(int bufferIndex, final double[] partial) {
            assert(partial.length == dimPartial);
            assert(partials != null);

            System.arraycopy(partial, 0, partials, dimPartial * bufferIndex, dimPartial);
        }

        @Override
        public void getPartial(int bufferIndex, final double[] partial) {
            assert(partial.length == dimPartial);
            assert(partials != null);

            System.arraycopy(partials, dimPartial * bufferIndex, partial, 0, dimPartial);
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
                final double remainder = remainders[rootBufferIndex * dimTrait + trait];

                logLikelihoods[trait] = logLike + remainder;

                if (incrementOuterProducts) {
                    int opo = dimTrait * dimTrait * trait;;

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
                    System.err.println("\t" + logLike + " " + (logLike + remainder));
//                    double[][] prec = new double[dimTrait][dimTrait];
//                    double[] mean = new double[dimTrait];
//                    for (int i = 0; i < dimTrait; ++i) {
//                        mean[i] = partials[rootOffset + i];
//                        for (int j = 0; j < dimTrait; ++j) {
//                            System.err.print(" " + diffusions[precisionOffset + i * dimTrait + j]);
//                            prec[i][j] = rootScalar * diffusions[precisionOffset + i * dimTrait + j];
//                        }
//                    }
//                    System.err.println("");
//                    MultivariateNormalDistribution mvn = new MultivariateNormalDistribution(new double[dimTrait],
//                            prec);
//                    System.err.println(mvn.logPdf(mean));
//                    System.err.println(mvn.getLogDet() + " ?= " + (dimTrait * Math.log(rootScalar) + precisionLogDet));
                }

                rootOffset += dimPartialForTrait;
                priorOffset += dimPartialForTrait;
            }
        }

        @Override
        public void updatePartials(final int[] operations, int operationCount, boolean incrementOuterProducts) {

            if (DEBUG) {
                System.err.println("Operations:");
            }

            int offset = 0;
            for (int op = 0; op < operationCount; ++op) {

                if (DEBUG) {
                    System.err.println("\t" + getOperationString(operations, offset));
                }

                updatePartial(
                        operations[offset + 0],
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
        public void updateDiffusionMatrices(int precisionIndex, final int[] probabilityIndices,
                                            final double[] edgeLengths, int updateCount) {

            if (DEBUG) {
                System.err.println("Matrices:");
            }

            for (int up = 0; up < updateCount; ++up) {

                if (DEBUG) {
                    System.err.println("\t" + probabilityIndices[up] + " <- " + edgeLengths[up]);
                }

                // TODO Currently only writtern for SCALAR model
                variances[dimMatrix * probabilityIndices[up]] = edgeLengths[up];
            }

            precisionOffset = dimTrait * dimTrait * precisionIndex;
            precisionLogDet = determinants[precisionIndex];
        }

        @Override
        public InstanceDetails getDetails() {
            return details;
        }

        // Internal storage
        private double[] partials;
        private double[] variances;
        private double[] remainders;
        private double[] diffusions;
        private double[] determinants;
        private int[] degreesOfFreedom;
        private double[] outerProducts;

        // Set during updateDiffusionMatrices() and used in updatePartials()
        private int precisionOffset;
        private double precisionLogDet;

        private static final boolean INLINE = true;

        private void updatePartial(
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
            final double vi = variances[imo];
            final double vj = variances[jmo];

            if (DEBUG) {
                System.err.println("i:");
                System.err.println("\tvar : " + variances[imo]);
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

                // Computer remainder at node k
                double remainder = 0.0;
                if (pi != 0 && pj != 0) {
                    if (INLINE) {

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

                    } else {
                        // incrementRemainderDensities(); // TODO
                    }

                    if (DEBUG) {
                        System.err.println("\t\tremainder: " + remainder);
                    }

                    if (incrementOuterProducts) {
                        int opo = dimTrait * dimTrait * trait;

                        final double remainderPrecision = pip * pjp / (pip + pjp);

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

                        degreesOfFreedom[trait] += 1; // incremenent degrees-of-freedom
                    }
                } // End if remainder

                // Accumulate remainder up tree and store
                remainders[kBuffer * dimTrait + trait] = remainder
                        + remainders[iBuffer * dimTrait + trait] + remainders[jBuffer * dimTrait + trait];

                if (DEBUG) {
                    System.err.println("\ttrait: " + trait);
                    System.err.println("\t\tprec: " + pi);
                    System.err.print("\t\tmean:");
                    for (int e = 0; e < dimTrait; ++e) {
                        System.err.print(" " + partials[ibo + e]);
                    }
                    System.err.println("");
                    System.err.println("");
                }

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
            variances = new double[dimMatrix * bufferCount];
            remainders = new double[bufferCount * dimTrait];

            diffusions = new double[dimTrait * dimTrait * diffusionCount];
            determinants = new double[diffusionCount];

            degreesOfFreedom = new int[numTraits];
            outerProducts = new double[dimTrait * dimTrait * numTraits];
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

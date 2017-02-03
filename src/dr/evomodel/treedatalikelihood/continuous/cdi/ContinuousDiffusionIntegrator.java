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
import org.ejml.alg.dense.decomposition.lu.LUDecompositionAlt_D64;
import org.ejml.alg.dense.linsol.lu.LinearSolverLu_D64;
import org.ejml.alg.dense.misc.UnrolledDeterminantFromMinor;
import org.ejml.alg.dense.misc.UnrolledInverseFromMinor;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

import java.util.Arrays;

import static dr.evomodel.treedatalikelihood.continuous.cdi.ContinuousDiffusionIntegrator.Multivariate.InversionResult.Code.FULLY_OBSERVED;
import static dr.evomodel.treedatalikelihood.continuous.cdi.ContinuousDiffusionIntegrator.Multivariate.InversionResult.Code.NOT_OBSERVED;
import static dr.evomodel.treedatalikelihood.continuous.cdi.ContinuousDiffusionIntegrator.Multivariate.InversionResult.Code.PARTIALLY_OBSERVED;
import static dr.evomodel.treedatalikelihood.continuous.cdi.ContinuousDiffusionIntegrator.Multivariate.wrap;

/**
 * @author Marc A. Suchard
 */
public interface ContinuousDiffusionIntegrator {
    int OPERATION_TUPLE_SIZE = 5;
    int NONE = -1;

    double LOG_SQRT_2_PI = 0.5 * Math.log(2 * Math.PI);

    void finalize() throws Throwable;

    void setPostOrderPartial(int bufferIndex, final double[] partial);

    void getPostOrderPartial(int bufferIndex, final double[] partial);

    void setPreOrderPartial(int bufferIndex, final double[] partial);

    void getPreOrderPartial(int bufferIndex, final double[] partial);

    void setWishartStatistics(final int[] degreesOfFreedom, final double[] outerProducts);

    void getWishartStatistics(final int[] degreesOfFreedom, final double[] outerProducts);

    void setDiffusionPrecision(int diffusionIndex, final double[] matrix, double logDeterminant);

    void updatePostOrderPartials(final int[] operations, int operationCount, boolean incrementOuterProducts);

    void updatePreOrderPartials(final int[] operations, int operationCount);

    InstanceDetails getDetails();

    void updateDiffusionMatrices(int precisionIndex, final int[] probabilityIndices, final double[] edgeLengths,
                                 int updateCount);

    void calculateRootLogLikelihood(int rootBufferIndex, int priorBufferIndex, double[] logLike,
                                    boolean incrementOuterProducts);

    boolean requireDataAugmentationForOuterProducts();

    // TODO Only send a list of operations
    void updatePreOrderPartial(int kp, int ip, int im, int jp, int jm);

    class Basic implements ContinuousDiffusionIntegrator {

        private int instance = -1;
        private InstanceDetails details = new InstanceDetails();

        private final PrecisionType precisionType;
        protected final int numTraits;
        protected final int dimTrait;
        protected final int bufferCount;
        protected final int diffusionCount;

        protected final int dimMatrix;
        protected final int dimPartialForTrait;
        protected final int dimPartial;

//        public Basic(final Basic base) {
//            this(base.precisionType, base.numTraits, base.dimTrait, base.bufferCount, base.diffusionCount);
//        }

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
        public boolean requireDataAugmentationForOuterProducts() {
            return false;
        }

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
            assert(partial.length == dimPartial);
            assert(partials != null);

            System.arraycopy(partials, dimPartial * bufferIndex, partial, 0, dimPartial);
        }

        @Override
        public void setPreOrderPartial(int bufferIndex, final double[] partial) {

            if (partial.length != dimPartial) {
                System.err.println("pl = " + partial.length);
                System.err.println("dp = " + dimPartial);
            }

            assert(partial.length == dimPartial);
            assert(prePartials != null);

            System.arraycopy(partial, 0, prePartials, dimPartial * bufferIndex, dimPartial);
        }

        @Override
        public void getPreOrderPartial(int bufferIndex, final double[] partial) {
            assert(partial.length == dimPartial);
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
//                System.exit(-1);
            }
        }

        @Override
        public void updatePreOrderPartials(final int[] operations, int operationCount) {
            throw new RuntimeException("Not yet implemented");
        }

        @Override
        public void updatePostOrderPartials(final int[] operations, int operationCount, boolean incrementOuterProducts) {

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
        protected double[] partials;
        protected double[] variances;
        protected double[] remainders;
        protected double[] diffusions;
        protected double[] determinants;
        protected int[] degreesOfFreedom;
        protected double[] outerProducts;
        protected double[] prePartials;

        // Set during updateDiffusionMatrices() and used in updatePartials()
        protected int precisionOffset;
        protected double precisionLogDet;

        protected static final boolean INLINE = true;

        @Override
        public void updatePreOrderPartial(
                final int kBuffer, // parent
                final int iBuffer, // node
                final int iMatrix,
                final int jBuffer, // sibling
                final int jMatrix) {
            throw new RuntimeException("Not yet implemented");
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

                        if (DEBUG && incrementOuterProducts) {
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
            variances = new double[dimMatrix * bufferCount];
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

    class Multivariate extends Basic {

        private static boolean DEBUG = false;

        public Multivariate(PrecisionType precisionType, int numTraits, int dimTrait, int bufferCount,
                            int diffusionCount) {
            super(precisionType, numTraits, dimTrait, bufferCount, diffusionCount);

            assert precisionType == PrecisionType.FULL;

            allocateStorage();
        }

        private void allocateStorage() {
            inverseDiffusions = new double[dimTrait * dimTrait * diffusionCount];
        }

        @Override
        public void setDiffusionPrecision(int precisionIndex, final double[] matrix, double logDeterminant) {
            super.setDiffusionPrecision(precisionIndex, matrix, logDeterminant);

            assert (inverseDiffusions != null);

            final int offset = dimTrait * dimTrait * precisionIndex;
            DenseMatrix64F precision = wrap(diffusions, offset, dimTrait, dimTrait);
            DenseMatrix64F variance = new DenseMatrix64F(dimTrait, dimTrait);
            CommonOps.invert(precision, variance);
            unwrap(variance, inverseDiffusions, offset);

            if (DEBUG) {
                System.err.println("At precision index: " + precisionIndex);
                System.err.println("precision: " + precision);
                System.err.println("variance : " + variance);
            }
        }

        private static void invert(final double[] source, final int sourceOffset,
                                   final double[] destination, final int destinationOffset,
                                   final int dim) {
//            throw new RuntimeException("Not yet implemented");
        }

        public static DenseMatrix64F wrap(final double[] source, final int offset,
                                          final int numRows, final int numCols) {
            double[] buffer = new double[numRows * numCols];
            return wrap(source, offset, numRows, numCols, buffer);
        }

        public static DenseMatrix64F wrap(final double[] source, final int offset,
                                                  final int numRows, final int numCols,
                                                  final double[] buffer) {
            System.arraycopy(source, offset, buffer, 0, numRows * numCols);
            return DenseMatrix64F.wrap(numRows, numCols, buffer);
        }


        public static void gatherRowsAndColumns(final DenseMatrix64F source, final DenseMatrix64F destination,
                                                          final int[] rowIndices, final int[] colIndices) {

            final int rowLength = rowIndices.length;
            final int colLength = colIndices.length;
            final double[] out = destination.getData();

            int index = 0;
            for (int i = 0; i < rowLength; ++i) {
                final int rowIndex = rowIndices[i];
                for (int j = 0; j < colLength; ++j) {
                    out[index] = source.unsafe_get(rowIndex, colIndices[j]);
                    ++index;
                }
            }
        }

        public static void scatterRowsAndColumns(final DenseMatrix64F source, final DenseMatrix64F destination,
                                                 final int[] rowIdices, final int[] colIndices, final boolean clear) {
            if (clear) {
                Arrays.fill(destination.getData(), 0.0);
            }

            final int rowLength = rowIdices.length;
            final int colLength = colIndices.length;
            final double[] in = source.getData();

            int index = 0;
            for (int i = 0; i < rowLength; ++i) {
                final int rowIndex = rowIdices[i];
                for (int j = 0; j < colLength; ++j) {
                    destination.unsafe_set(rowIndex, colIndices[j], in[index]);
                    ++index;
                }
            }
        }


        public static void unwrap(final DenseMatrix64F source, final double[] destination, final int offset) {
            System.arraycopy(source.getData(), 0, destination, offset, source.getNumElements());
        }

        public static boolean anyDiagonalInfinities(DenseMatrix64F source) {
            boolean anyInfinities = false;
            for (int i = 0; i < source.getNumCols() && !anyInfinities; ++i) {
                if (Double.isInfinite(source.unsafe_get(i, i))) {
                    anyInfinities = true;
                }
            }
            return anyInfinities;
        }

        public static boolean allFiniteDiagonals(DenseMatrix64F source) {
            boolean allFinite = true;

            final int length = source.getNumCols();
            for (int i = 0; i < length; ++i) {
                allFinite &= !Double.isInfinite(source.unsafe_get(i, i));
            }
            return allFinite;
        }

        public static int countFiniteDiagonals(DenseMatrix64F source) {
            final int length = source.getNumCols();

            int count = 0;
            for (int i = 0; i < length; ++i) {
                final double d = source.unsafe_get(i, i);
                if (!Double.isInfinite(d)) {
                    ++count;
                }
            }
            return count;
        }

        public static int countZeroDiagonals(DenseMatrix64F source) {
            final int length = source.getNumCols();

            int count = 0;
            for (int i = 0; i < length; ++i) {
                final double d = source.unsafe_get(i, i);
                if (d == 0.0) {
                    ++count;
                }
            }
            return count;
        }

        public static void getFiniteDiagonalIndices(final DenseMatrix64F source, final int[] indices) {
            final int length = source.getNumCols();

            int index = 0;
            for (int i = 0; i < length; ++i) {
                final double d = source.unsafe_get(i, i);
                if (!Double.isInfinite(d)) {
                    indices[index] = i;
                    ++index;
                }
            }
        }

        public static int countFiniteNonZeroDiagonals(DenseMatrix64F source) {
            final int length = source.getNumCols();

            int count = 0;
            for (int i = 0; i < length; ++i) {
                final double d = source.unsafe_get(i, i);
                if (!Double.isInfinite(d) && d != 0.0) {
                    ++count;
                }
            }
            return count;
        }

        public static void getFiniteNonZeroDiagonalIndices(final DenseMatrix64F source, final int[] indices) {
            final int length = source.getNumCols();

            int index = 0;
            for (int i = 0; i < length; ++i) {
                final double d = source.unsafe_get(i, i);
                if (!Double.isInfinite(d) && d != 0.0) {
                    indices[index] = i;
                    ++index;
                }
            }
        }

//        public static void invertWithInfiniteDiagonals(DenseMatrix64F source, DenseMatrix64F destination) {
//
//            if (anyDiagonalInfinities(source)) {
//                throw new RuntimeException("Not yet implemented");
//            } else {
//                CommonOps.invert(source, destination);
//            }
//        }
//
//        public static DenseMatrix64F invertWithInfiniteDiagonals(DenseMatrix64F source) {
//            DenseMatrix64F destination = new DenseMatrix64F(source.getNumRows(), source.getNumCols());
//            invertWithInfiniteDiagonals(source, destination);
//            return destination;
//        }
//
//        public static void addToDiagonal(DenseMatrix64F source, DenseMatrix64F destination, double increment) {
//            System.arraycopy(source.getData(), 0, destination.getData(), 0, source.getNumElements());
//            addToDiagonal(destination, increment);
//        }

        public static void addToDiagonal(DenseMatrix64F source, double increment) {
            final int width = source.getNumRows();
            for (int i = 0; i < width; ++i) {
                source.unsafe_set(i,i, source.unsafe_get(i, i) + increment);
            }
        }



        public static class InversionResult {

            enum Code {
                FULLY_OBSERVED,
                NOT_OBSERVED,
                PARTIALLY_OBSERVED
            }

            InversionResult(Code code, int dim, double determinant) {
                this.code = code;
                this.dim = dim;
                this.determinant = determinant;
            }

            final public Code getReturnCode() {
                return code;
            }

            final public int getEffectiveDimension() {
                return dim;
            }

            final public double getDeterminant() {
                return determinant;
            }

            public String toString() {
                return code + ":" + dim;
            }

            final private Code code;
            final private int dim;
            final private double determinant;
        }

        public static double det(DenseMatrix64F mat) {
            int numCol = mat.getNumCols();
            int numRow = mat.getNumRows();
            if(numCol != numRow) {
                throw new IllegalArgumentException("Must be a square matrix.");
            } else if(numCol <= 6) {
                return numCol >= 2? UnrolledDeterminantFromMinor.det(mat):mat.get(0);
            } else {
                LUDecompositionAlt_D64 alg = new LUDecompositionAlt_D64();
                if(alg.inputModified()) {
                    mat = mat.copy();
                }

                return !alg.decompose(mat)?0.0D:alg.computeDeterminant().real;
            }
        }

        public static double invertAndGetDeterminant(DenseMatrix64F mat, DenseMatrix64F result) {

            final int numCol = mat.getNumCols();
            final int numRow = mat.getNumRows();
            if (numCol != numRow) {
                throw new IllegalArgumentException("Must be a square matrix.");
            }

            if (numCol <= 5) {

                if (numCol >= 2) {
                    UnrolledInverseFromMinor.inv(mat, result);
                } else {
                    result.set(0, 1.0D / mat.get(0));
                }

                return numCol >= 2 ?
                        UnrolledDeterminantFromMinor.det(mat) :
                        mat.get(0);
            } else {

                LUDecompositionAlt_D64 alg = new LUDecompositionAlt_D64();
                LinearSolverLu_D64 solver = new LinearSolverLu_D64(alg);
                if (solver.modifiesA()) {
                    mat = mat.copy();
                }

                if (!solver.setA(mat)) {
                    return Double.NaN;
                }

                solver.invert(result);

                return alg.computeDeterminant().real;

            }
        }


        public static InversionResult safeInvert(DenseMatrix64F source, DenseMatrix64F destination, boolean getDeterminant) {

            final int dim = source.getNumCols();
            final int finiteCount = countFiniteNonZeroDiagonals(source);
            double det = 0;

            if (finiteCount == dim) {
                if (getDeterminant) {
                    det = invertAndGetDeterminant(source, destination);
                } else {
                    CommonOps.invert(source, destination);
                }
                return new InversionResult(FULLY_OBSERVED, dim, det);
            } else {
                if (finiteCount == 0) {
                    Arrays.fill(destination.getData(), 0);
                    return new InversionResult(NOT_OBSERVED, 0, 0);
                } else {
                    final int[] finiteIndices = new int[finiteCount];
                    getFiniteNonZeroDiagonalIndices(source, finiteIndices);

                    final DenseMatrix64F subSource = new DenseMatrix64F(finiteCount, finiteCount);
                    gatherRowsAndColumns(source, subSource, finiteIndices, finiteIndices);

                    final DenseMatrix64F inverseSubSource = new DenseMatrix64F(finiteCount, finiteCount);
                    if (getDeterminant) {
                        det = invertAndGetDeterminant(subSource, inverseSubSource);
                    } else {
                        CommonOps.invert(subSource, inverseSubSource);
                    }

                    scatterRowsAndColumns(inverseSubSource, destination, finiteIndices, finiteIndices, true);

                    return new InversionResult(PARTIALLY_OBSERVED, finiteCount, det);
                }
            }
        }

        public static void weightedAverage(final WrappedVector mi,
                                           final DenseMatrix64F Pi,
                                           final WrappedVector mj,
                                           final DenseMatrix64F Pj,
                                           final WrappedVector mk,
                                           final DenseMatrix64F Vk,
                                           final int dimTrait) {
            final double[] tmp = new double[dimTrait];
            for (int g = 0; g < dimTrait; ++g) {
                double sum = 0.0;
                for (int h = 0; h < dimTrait; ++h) {
                    sum += Pi.unsafe_get(g, h) * mi.get(h);
                    sum += Pj.unsafe_get(g, h) * mj.get(h);
                }
                tmp[g] = sum;
            }
            for (int g = 0; g < dimTrait; ++g) {
                double sum = 0.0;
                for (int h = 0; h < dimTrait; ++h) {
                    sum += Vk.unsafe_get(g, h) * tmp[h];
                }
                mk.set(g, sum);
            }
        }

        public static void weightedAverage(final double[] ipartial,
                                           final int ibo,
                                           final DenseMatrix64F Pi,
                                           final double[] jpartial,
                                           final int jbo,
                                           final DenseMatrix64F Pj,
                                           final double[] kpartial,
                                           final int kbo,
                                           final DenseMatrix64F Vk,
                                           final int dimTrait) {
            final double[] tmp = new double[dimTrait];
            for (int g = 0; g < dimTrait; ++g) {
                double sum = 0.0;
                for (int h = 0; h < dimTrait; ++h) {
                    sum += Pi.unsafe_get(g, h) * ipartial[ibo + h];
                    sum += Pj.unsafe_get(g, h) * jpartial[jbo + h];
                }
                tmp[g] = sum;
            }
            for (int g = 0; g < dimTrait; ++g) {
                double sum = 0.0;
                for (int h = 0; h < dimTrait; ++h) {
                    sum += Vk.unsafe_get(g, h) * tmp[h];
                }
                kpartial[kbo + g] = sum;
            }
        }

        @Override
        public boolean requireDataAugmentationForOuterProducts() {
            return true;
        }

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
            final double vi = variances[imo];
            final double vj = variances[jmo];

            final DenseMatrix64F Vd = wrap(inverseDiffusions, precisionOffset, dimTrait, dimTrait);

            if (DEBUG) {
                System.err.println("updatePreOrderPartial for node " + iBuffer);
//                System.err.println("variance diffusion: " + Vd);
                System.err.println("\tvi: " + vi + " vj: " + vj);
//                System.err.println("precisionOffset = " + precisionOffset);
            }

            // For each trait // TODO in parallel
            for (int trait = 0; trait < numTraits; ++trait) {

                // A. Get current precision of k and j
                final DenseMatrix64F Pk = wrap(prePartials, kbo + dimTrait, dimTrait, dimTrait);
//                final DenseMatrix64F Pj = wrap(partials, jbo + dimTrait, dimTrait, dimTrait);

//                final DenseMatrix64F Vk = wrap(prePartials, kbo + dimTrait + dimTrait * dimTrait, dimTrait, dimTrait);
                final DenseMatrix64F Vj = wrap(partials, jbo + dimTrait + dimTrait * dimTrait, dimTrait, dimTrait);

                // B. Inflate variance along sibling branch using matrix inversion
                final DenseMatrix64F Vjp = new DenseMatrix64F(dimTrait, dimTrait);
                CommonOps.add(Vj, vj, Vd, Vjp);

                final DenseMatrix64F Pjp = new DenseMatrix64F(dimTrait, dimTrait);
                InversionResult cj = safeInvert(Vjp, Pjp, false);

                final DenseMatrix64F Pip = new DenseMatrix64F(dimTrait, dimTrait);
                CommonOps.add(Pk, Pjp, Pip);

                final DenseMatrix64F Vip = new DenseMatrix64F(dimTrait, dimTrait);
                InversionResult cip = safeInvert(Pip, Vip, false);

                // C. Compute prePartial mean
                final double[] tmp = new double[dimTrait];
                for (int g = 0; g < dimTrait; ++g) {
                    double sum = 0.0;
                    for (int h = 0; h < dimTrait; ++h) {
                        sum += Pk.unsafe_get(g, h) * prePartials[kbo + h]; // Read parent
                        sum += Pjp.unsafe_get(g, h) * partials[jbo + h];   // Read sibling
                    }
                    tmp[g] = sum;
                }
                for (int g = 0; g < dimTrait; ++g) {
                    double sum = 0.0;
                    for (int h = 0; h < dimTrait; ++h) {
                        sum += Vip.unsafe_get(g, h) * tmp[h];
                    }
                    prePartials[ibo + g] = sum; // Write node
                }

                // C. Inflate variance along node branch
                final DenseMatrix64F Vi = Vip;
                CommonOps.add(vi, Vd, Vip, Vi);

                final DenseMatrix64F Pi = new DenseMatrix64F(dimTrait, dimTrait);
                InversionResult ci = safeInvert(Vi, Pi, false);

                // X. Store precision results for node
                unwrap(Pi, prePartials, ibo + dimTrait);
                unwrap(Vi, prePartials, ibo + dimTrait + dimTrait * dimTrait);

                if (DEBUG) {
                    System.err.println("trait: " + trait);
                    System.err.println("pM: " + new WrappedVector.Raw(prePartials, kbo, dimTrait));
                    System.err.println("pP: " + Pk);
                    System.err.println("sM: " + new WrappedVector.Raw(partials, jbo, dimTrait));
                    System.err.println("sV: " + Vj);
                    System.err.println("sVp: " + Vjp);
                    System.err.println("sPp: " + Pjp);
                    System.err.println("Pip: " + Pip);
                    System.err.println("cM: " + new WrappedVector.Raw(prePartials, ibo, dimTrait));
                    System.err.println("cV: " + Vi);
                }

                // Get ready for next trait
                kbo += dimPartialForTrait;
                ibo += dimPartialForTrait;
                jbo += dimPartialForTrait;
            }
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

            final DenseMatrix64F Vd = wrap(inverseDiffusions, precisionOffset, dimTrait, dimTrait);

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

                // Increase variance along the branches i -> k and j -> k

                // A. Get current precision of i and j
                final double lpi = partials[ibo + dimTrait + 2 * dimTrait * dimTrait];
                final double lpj = partials[jbo + dimTrait + 2 * dimTrait * dimTrait];

                final DenseMatrix64F Pi = wrap(partials, ibo + dimTrait, dimTrait, dimTrait);
                final DenseMatrix64F Pj = wrap(partials, jbo + dimTrait, dimTrait, dimTrait);

                final DenseMatrix64F Vi = wrap(partials, ibo + dimTrait + dimTrait * dimTrait, dimTrait, dimTrait);
                final DenseMatrix64F Vj = wrap(partials, jbo + dimTrait + dimTrait * dimTrait, dimTrait, dimTrait);

                // B. Integrate along branch using two matrix inversions
                final double lpip = Double.isInfinite(lpi) ?
                        1.0 / vi : lpi / (1.0 + lpi * vi);
                final double lpjp = Double.isInfinite(lpj) ?
                        1.0 / vj : lpj / (1.0 + lpj * vj);

                final DenseMatrix64F Vip = new DenseMatrix64F(dimTrait, dimTrait);
                final DenseMatrix64F Vjp = new DenseMatrix64F(dimTrait, dimTrait);

                CommonOps.add(Vi, vi, Vd, Vip);
                CommonOps.add(Vj, vj, Vd, Vjp);

                final DenseMatrix64F Pip = new DenseMatrix64F(dimTrait, dimTrait);
                final DenseMatrix64F Pjp = new DenseMatrix64F(dimTrait, dimTrait);

                InversionResult ci = safeInvert(Vip, Pip, true);
                InversionResult cj = safeInvert(Vjp, Pjp, true);

                // Compute partial mean and precision at node k

                // A. Partial precision and variance (for later use) using one matrix inversion
                final double lpk = lpip + lpjp;

                final DenseMatrix64F Pk = new DenseMatrix64F(dimTrait, dimTrait);

                CommonOps.add(Pip, Pjp, Pk);

                final DenseMatrix64F Vk = new DenseMatrix64F(dimTrait, dimTrait);
                InversionResult ck = safeInvert(Pk, Vk, true);

                // B. Partial mean
//                for (int g = 0; g < dimTrait; ++g) {
//                    partials[kbo + g] = (pip * partials[ibo + g] + pjp * partials[jbo + g]) / pk;
//                }

                final double[] tmp = new double[dimTrait];
                for (int g = 0; g < dimTrait; ++g) {
                    double sum = 0.0;
                    for (int h = 0; h < dimTrait; ++h) {
                        sum += Pip.unsafe_get(g, h) * partials[ibo + h];
                        sum += Pjp.unsafe_get(g, h) * partials[jbo + h];
                    }
                    tmp[g] = sum;
                }
                for (int g = 0; g < dimTrait; ++g) {
                    double sum = 0.0;
                    for (int h = 0; h < dimTrait; ++h) {
                        sum += Vk.unsafe_get(g, h) * tmp[h];
                    }
                    partials[kbo + g] = sum;
                }

                // C. Store precision
                partials[kbo + dimTrait + 2 * dimTrait * dimTrait] = lpk;

                unwrap(Pk, partials, kbo + dimTrait);
                unwrap(Vk, partials, kbo + dimTrait + dimTrait * dimTrait);

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
                    System.err.println("Vip: " + Vip);
                    System.err.println("Pjp: " + Pjp);
                    System.err.println("Vjp: " + Vjp);
                }

                if (!(ci.getReturnCode() == NOT_OBSERVED || cj.getReturnCode() == NOT_OBSERVED)) {
//                if (ci == InversionReturnCode.FULLY_OBSERVED && cj == InversionReturnCode.FULLY_OBSERVED) {
                    // TODO Fix for partially observed
//                if (pi != 0 && pj != 0) {
//

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

                    final DenseMatrix64F Vt = new DenseMatrix64F(dimTrait, dimTrait);
                    CommonOps.add(Vip, Vjp, Vt);

                    if (DEBUG) {
                        System.err.println("Vt: " + Vt);
                    }

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
                        System.err.println("\t\tremainder: " + remainder);
//                        System.exit(-1);
                    }

                    if (incrementOuterProducts) {

                        assert (false);

                        final DenseMatrix64F Pt = new DenseMatrix64F(dimTrait, dimTrait);
                        InversionResult ct = safeInvert(Vt, Pt, false);

                        int opo = dimTrait * dimTrait * trait;
                        int opd = precisionOffset;

                        for (int g = 0; g < dimTrait; ++g) {
                            final double ig = partials[ibo + g];
                            final double jg = partials[jbo + g];

                            for (int h = 0; h < dimTrait; ++h) {
                                final double ih = partials[ibo + h];
                                final double jh = partials[jbo + h];

                                outerProducts[opo] += (ig - jg) * (ih - jh)
//                                        * Pt.unsafe_get(g, h)
//                                        * Pk.unsafe_get(g, h)

//                                        / diffusions[opd];
                                 // * pip * pjp / (pip + pjp);
                                        * lpip * lpjp / (lpip + lpjp);
                                ++opo;
                                ++opd;
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

            final DenseMatrix64F Vd = wrap(inverseDiffusions, precisionOffset, dimTrait, dimTrait);

            // TODO For each trait in parallel
            for (int trait = 0; trait < numTraits; ++trait) {

                final DenseMatrix64F Proot = wrap(partials, rootOffset + dimTrait, dimTrait, dimTrait);
                final DenseMatrix64F Pprior = wrap(partials, priorOffset + dimTrait, dimTrait, dimTrait);

                final DenseMatrix64F Vroot = wrap(partials, rootOffset + dimTrait + dimTrait * dimTrait, dimTrait, dimTrait);
                final DenseMatrix64F Vprior = wrap(partials, priorOffset + dimTrait + dimTrait * dimTrait, dimTrait, dimTrait);

                // TODO Block below is for the conjugate prior ONLY
                {
                    final DenseMatrix64F Vtmp = new DenseMatrix64F(dimTrait, dimTrait);
                    CommonOps.mult(Vd, Vprior, Vtmp);
                    Vprior.set(Vtmp);
                }

                final DenseMatrix64F Vtotal = new DenseMatrix64F(dimTrait, dimTrait);
                CommonOps.add(Vroot, Vprior, Vtotal);

                final DenseMatrix64F Ptotal = new DenseMatrix64F(dimTrait, dimTrait);
                CommonOps.invert(Vtotal, Ptotal);  // TODO Can return determinant at same time to avoid extra QR decomp

                double SS = 0;
                for (int g = 0; g < dimTrait; ++g) {
                    final double gDifference = partials[rootOffset + g] - partials[priorOffset + g];

                    for (int h = 0; h < dimTrait; ++h) {
                        final double hDifference = partials[rootOffset + h] - partials[priorOffset + h];

                        SS += gDifference * Ptotal.unsafe_get(g, h) * hDifference;
                    }
                }

                final double logLike = -dimTrait * LOG_SQRT_2_PI - 0.5 * Math.log(CommonOps.det(Vtotal)) - 0.5 * SS;

                final double remainder = remainders[rootBufferIndex * numTraits + trait];
                logLikelihoods[trait] = logLike + remainder;

                if (incrementOuterProducts) {
                    int opo = dimTrait * dimTrait * trait;
                    int opd = precisionOffset;

                    double rootScalar = partials[rootOffset + dimTrait + 2 * dimTrait * dimTrait];
                    final double priorScalar = partials[priorOffset + dimTrait];

                    if (!Double.isInfinite(priorScalar)) {
                        rootScalar = rootScalar * priorScalar / (rootScalar + priorScalar);
                    }

                    for (int g = 0; g < dimTrait; ++g) {
                        final double gDifference = partials[rootOffset + g] - partials[priorOffset + g];

                        for (int h = 0; h < dimTrait; ++h) {
                            final double hDifference = partials[rootOffset + h] - partials[priorOffset + h];

                            outerProducts[opo] += gDifference * hDifference
//                                    * Ptotal.unsafe_get(g, h) / diffusions[opd];
                            * rootScalar;
                            ++opo;
                            ++opd;
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
                    System.err.println("Proot: " + Proot);
                    System.err.println("Vroot: " + Vroot);
                    System.err.println("Pprior: " + Pprior);
                    System.err.println("Vprior: " + Vprior);
                    System.err.println("Ptotal: " + Ptotal);
                    System.err.println("Vtotal: " + Vtotal);
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
//                System.exit(-1);
            }
        }

        private double[] inverseDiffusions;
    }
}

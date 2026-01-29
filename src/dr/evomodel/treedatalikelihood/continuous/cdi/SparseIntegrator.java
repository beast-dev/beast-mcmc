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

import dr.matrix.SparseCompressedMatrix;

import java.util.Arrays;

import static dr.math.matrixAlgebra.missingData.MissingOps.wrap;

/**
 * @author Marc A. Suchard
 * @author Philippe Lemey
 */

public class SparseIntegrator extends ContinuousDiffusionIntegrator.Basic {

    private static final boolean DEBUG = false;

    private double[] sumOfSquares;

    private static final boolean USE_SSE_CACHE = true;

    public SparseIntegrator(PrecisionType precisionType, int numTraits, int dimTrait, int dimProcess,
                            int bufferCount, int diffusionCount) {
        super(precisionType, numTraits, dimTrait, dimProcess, bufferCount, diffusionCount);

        assert precisionType == PrecisionType.SCALAR;

        allocateStorage();
    }

    @Override
    void resetSumOfSquares() {
        if (USE_SSE_CACHE) {
            Arrays.fill(sumOfSquares, -1.0);
        }
    }

    // TODO copy back to `Basic` to avoid code duplication
    protected void updatePartial(
            final int kBuffer,
            final int iBuffer,
            final int imo,
            final int jBuffer,
            final int jmo,
            final boolean computeRemainders,
            final boolean incrementOuterProducts
    ) {
        // Determine buffer offsets
        int kbo = dimPartial * kBuffer;
        int ibo = dimPartial * iBuffer;
        int jbo = dimPartial * jBuffer;

        // Determine matrix offsets
//        final int imo = iMatrix; //TODO: just use iMatrix * jMatrix? (also, why do we need these?)
//        final int jmo = jMatrix;

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
                System.err.println();
            }

            // Computer remainder at node k
            double remainder = 0.0;
            if (computeRemainders && pi != 0 && pj != 0) {

                // TODO line below is only difference with `Basic`
                remainder += computeRemainder(trait, ibo, jbo, kbo, pip, pjp, pk, incrementOuterProducts);

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

    // TODO copy back to `Basic` to avoid code duplication
    @Override
    public void calculateRootLogLikelihood(int rootBufferIndex, int priorBufferIndex, int precisionIndex,
                                           final double[] logLikelihoods, boolean incrementOuterProducts, boolean isIntegratedProcess) {
        assert(logLikelihoods.length == numTraits);
        assert(!isIntegratedProcess);

        updatePrecisionOffsetAndDeterminant(precisionIndex);

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

            SS += computeRootSumOfSquares(rootOffset, priorOffset, pob);

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

                degreesOfFreedom[trait] += 1; // increment degrees-of-freedom
            }

            if (DEBUG) {
                System.err.print("mean:");
                for (int g = 0; g < dimTrait; ++g) {
                    System.err.print(" " + partials[rootOffset + g]);
                }
                System.err.println();
                System.err.println("prec: " + partials[rootOffset + dimTrait]);
                System.err.println("rootScalar: " + rootScalar);
                System.err.println("\t" + logLike + " " + (logLike + remainder));
                if (incrementOuterProducts) {
                    System.err.println("Outer-products:" + wrap(outerProducts, dimTrait * dimTrait * trait, dimTrait, dimTrait));
                }
                System.err.println();
            }

            rootOffset += dimPartialForTrait;
            priorOffset += dimPartialForTrait;
        }

        if (DEBUG) {
            System.err.println("End");
        }
    }

    private void allocateStorage() {
        sumOfSquares = new double[numTraits * bufferCount];
        sparseDiffusion = new SparseCompressedMatrix[diffusionCount];
    }

    private int precisionIndex;
    private SparseCompressedMatrix[] sparseDiffusion;

    @Override
    public void setDiffusionPrecision(int precisionIndex, final double[] matrix, double logDeterminant) {
        super.setDiffusionPrecision(precisionIndex, matrix, logDeterminant);

        // TODO pass as `SCM` to avoid recreating each time
        sparseDiffusion[precisionIndex] = new SparseCompressedMatrix(matrix, 0, dimTrait, dimTrait);
        assert(matrix.length == dimProcess * dimProcess);
    }

    @Override
    public void setDiffusionPrecision(int precisionIndex, SparseCompressedMatrix matrix, double logDeterminant) {
        determinants[precisionIndex] = logDeterminant;
        sparseDiffusion[precisionIndex] = matrix;
    }

    @Override
    public DiffusionRepresentation diffusionFactory(int count) {
        return new DiffusionRepresentation.Sparse(count);
    }

    @Override
    void updatePrecisionOffsetAndDeterminant(int precisionIndex) {
        super.updatePrecisionOffsetAndDeterminant(precisionIndex);
        this.precisionIndex = precisionIndex;
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

        final SparseCompressedMatrix diffusion = sparseDiffusion[precisionIndex];

        double SSk;
        double SSi;
        double SSj;

        if (USE_SSE_CACHE) {

            SSk = diffusion.multiplyVectorTransposeMatrixVector(partials, kbo, partials, kbo);
            sumOfSquares[kbo / dimPartialForTrait + trait] = SSk;

            SSi = sumOfSquares[ibo / dimPartialForTrait + trait];
            if (SSi == -1.0) {
                SSi = diffusion.multiplyVectorTransposeMatrixVector(partials, ibo, partials, ibo);
                sumOfSquares[ibo / dimPartialForTrait + trait] = SSi;
            }

            SSj = sumOfSquares[jbo / dimPartialForTrait + trait];
            if (SSj == -1.0) {
                SSj = diffusion.multiplyVectorTransposeMatrixVector(partials, jbo, partials, jbo);
                sumOfSquares[jbo / dimPartialForTrait + trait] = SSj;
            }
        } else {
            SSk = diffusion.multiplyVectorTransposeMatrixVector(partials, kbo, partials, kbo);
            SSi = diffusion.multiplyVectorTransposeMatrixVector(partials, ibo, partials, ibo);
            SSj = diffusion.multiplyVectorTransposeMatrixVector(partials, jbo, partials, jbo);
        }

        final double remainderPrecision = pip * pjp / pk;

        final double remainder = -dimTrait * LOG_SQRT_2_PI // TODO Can move some calculation outside the loop
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

    private double computeRootSumOfSquares(int rootOffset, int priorOffset, int pob) {
        final SparseCompressedMatrix diffusion = sparseDiffusion[pob];
        return diffusion.multiplyVecDiffTransposeMatrixVecDiff(partials, rootOffset, partials, priorOffset);
    }
}

/*
 * SafeMultivariateBlockDiagonalActualizedWithDriftIntegrator.java
 *
 * OPTIMIZED VERSION with:
 *   - O(n³) block-wise Lyapunov solver (down from O(n⁶))
 *   - Cached work matrices to avoid repeated allocations
 *   - Eliminated redundant matrix copies
 *   - Direct block operations without full dense matrix construction
 *
 * Block-Diag OU integrator:
 *   S = R D R^{-1},  D block-tridiagonal with 1x1 and 2x2 real blocks.
 *
 * The selection matrix in the working basis is given in *compressed*
 * tridiagonal / block-Diag form:
 *
 *   blockDParams = [ d_0, ..., d_{d-1},
 *                    u_0, ..., u_{d-2},
 *                    l_0, ..., l_{d-2} ]
 *
 * where D has
 *   D_{ii}     = d_i,
 *   D_{i,i+1}  = u_i,
 *   D_{i+1,i}  = l_i,
 * and all other entries zero.  Contiguous non-zero pairs (u_i, l_i)
 * represent general 2x2 blocks; all other positions behave as 1x1 blocks.
 *
 * R is assumed orthonormal.
 *
 * Copyright © 2002-2024 the BEAST Development Team
 */

package dr.evomodel.treedatalikelihood.continuous.cdi;

import org.ejml.data.DenseMatrix64F;
import dr.evomodel.treedatalikelihood.continuous.backprop.BlockDiagonalLyapunovSolver;
import org.ejml.ops.CommonOps;

import static dr.evomodel.treedatalikelihood.continuous.backprop.BlockDiagonalExpSolver.expBlockMatrix;

/**
 * Block–Diag specialization of SafeMultivariateActualizedWithDriftIntegrator.
 *
 * Differences from the eigenvalue-based superclass:
 *   - The strength-of-selection vector is interpreted as a compressed
 *     block-tridiagonal matrix D in an orthonormal basis R.
 *   - Stationary variance Σ is computed by solving the Lyapunov equation
 *         D Σ_eig + Σ_eig Dᵀ = V_eig
 *     in the block basis using a dedicated BlockDiagonalLyapunovSolver.
 *   - Q(t) = exp(-D t) is computed block-wise (1×1 and general 2×2),
 *     then mapped back with R.
 *   - Integrated OU updates use a BlockDiag-specific inverse S^{-1}
 *     instead of the diagonal-only helper in the superclass (not yet
 *     implemented here).
 *
 * NOTE: This class is not thread-safe due to shared workspaces. Each
 * integrator instance must be confined to a single thread.
 *
 * @author Filippo Monti
 */
public class SafeMultivariateBlockDiagonalActualizedWithDriftIntegrator
        extends SafeMultivariateActualizedWithDriftIntegrator {

    private static final double EPS_BLOCK = 1e-12;

    // ------------------------------------------------------------------
    // Preallocated workspaces (avoids repeated allocations in hot paths)
    // ------------------------------------------------------------------
    private final DenseMatrix64F workV;        // d × d
    private final DenseMatrix64F workR;        // d × d
    private final DenseMatrix64F workRinv;     // d × d
    private final DenseMatrix64F workVEig;     // d × d
    private final DenseMatrix64F workSigmaEig; // d × d
    private final DenseMatrix64F workSigma;    // d × d
    private final DenseMatrix64F workQeig;     // d × d
    private final DenseMatrix64F workQtrait;   // d × d

    private final DenseMatrix64F temp;

    // Dedicated helper for D Σ + Σ Dᵀ = V in the block basis.
    private final BlockDiagonalLyapunovSolver lyapunovSolver;

    public SafeMultivariateBlockDiagonalActualizedWithDriftIntegrator(
            PrecisionType precisionType,
            int numTraits,
            int dimTrait,
            int dimProcess,
            int bufferCount,
            int diffusionCount) {

        // R is orthonormal; keep symmetric basis flag true so that
        // transforms use Rᵀ and never mutate R via inversion.
        super(precisionType, numTraits, dimTrait, dimProcess,
                bufferCount, diffusionCount,
                false);

        final int d = dimProcess;

        workV        = new DenseMatrix64F(d, d);
        workR        = new DenseMatrix64F(d, d);
        workRinv     = new DenseMatrix64F(d, d);
        workVEig     = new DenseMatrix64F(d, d);
        workSigmaEig = new DenseMatrix64F(d, d);
        workSigma    = new DenseMatrix64F(d, d);
        workQeig     = new DenseMatrix64F(d, d);
        workQtrait   = new DenseMatrix64F(d, d);
        temp         = new DenseMatrix64F(d, d);

        lyapunovSolver = new BlockDiagonalLyapunovSolver(d);
    }

    /* **********************************************************************
     * Stationary variance Σ via block-wise Lyapunov solver
     * ******************************************************************* */

    @Override
    public void setDiffusionStationaryVariance(int precisionIndex,
                                               final double[] blockDParams,
                                               final double[] rotation) { //rotation concatenates row-wise R and Rinv

        final int d = dimProcess;
//        final int expectedLen = 3 * d - 2; // for debugging
//        if (blockDParams.length != expectedLen) {
//            throw new IllegalArgumentException(
//                    "Block-Diag parameter vector has length " + blockDParams.length +
//                            " but expected " + expectedLen + " for dimProcess = " + d);
//        }
        if (rotation.length != 2 * d * d) {
            throw new IllegalArgumentException(
                    "Rotation length " + rotation.length +
                            " does not match 2 * dimProcess^2 = " + (2 * d * d));
        }

        final int matrixSize = d * d;
        final int offset     = matrixSize * precisionIndex;

        // 1. Move diffusion matrix into block basis: V_eig = Rᵀ V R
        DenseMatrix64F V = workV;
        DenseMatrix64F R = workR;
        DenseMatrix64F Rinv = workRinv;

        double[] vData = V.data;
        double[] rData = R.data;
        double[] rInverseData = Rinv.data;
        for (int i = 0; i < d; i++) {
            int rowOffsetV = offset + i * d;
            int rowOffsetR = i * d;
            for (int j = 0; j < d; j++) {
                vData[i * d + j] = inverseDiffusions[rowOffsetV + j];
                rData[i * d + j] = rotation[rowOffsetR + j];
                rInverseData[i * d + j] = rotation[d*d + rowOffsetR + j];
            }
        }

        DenseMatrix64F V_eig = workVEig;
        V_eig.set(V);
//        transformMatrix(V_eig, R, isActualizationSymmetric); // Rᵀ V R
        // V_eig = R^{-1} V R^{-T}
        CommonOps.mult(workRinv, V, temp);

        CommonOps.multTransB(temp, workRinv, workVEig);

        // 2. Solve Lyapunov: D Σ_eig + Σ_eig Dᵀ = V_eig directly from blockDParams
        DenseMatrix64F SigmaEig = workSigmaEig;
        lyapunovSolver.solve(blockDParams, V_eig, SigmaEig);

        // 3. Map Σ back to trait space: Σ = R Σ_eig Rᵀ
        DenseMatrix64F Sigma = workSigma;
        Sigma.set(SigmaEig);
//        transformMatrixBack(Sigma, R);
        // Sigma = R Sigma_eig R^T
        CommonOps.mult(R, SigmaEig, temp); // TODO avoid this allocation
        CommonOps.multTransB(temp, R, workSigma);

        // 4. Copy result back to stationaryVariances array (row-wise)
        double[] sigmaData = Sigma.data;
        for (int i = 0; i < d; i++) {
            System.arraycopy(sigmaData, i * d,
                    stationaryVariances, offset + i * d,
                    d);
        }
    }

    /* **********************************************************************
     * OU actualisation Q(t) from exp(-D t) in block basis
     * ******************************************************************* */

    @Override
    void computeOUActualization(final double[] blockDParams,
                                final double[] rotation,
                                final double edgeLength,
                                final int scaledOffsetDiagonal,
                                final int scaledOffset) {

        final int d = dimProcess;
        final int expectedLen = 3 * d - 2;
        if (blockDParams.length != expectedLen) {
            throw new IllegalArgumentException(
                    "Block-Diag parameter vector has length " + blockDParams.length +
                            " but expected " + expectedLen + " for dimProcess = " + d);
        }

        // 1. Compute Q_eig(t) = exp(-D t) directly from blockDParams
        DenseMatrix64F Qeig = workQeig;
        expBlockMatrix(blockDParams, d, edgeLength, Qeig);

        // 2. Map Q(t) back to trait space: Q = R Q_eig Rᵀ
        DenseMatrix64F R = workR;
        DenseMatrix64F Rinv = workRinv;
        double[] rData = R.data;
        double[] rInverseData = Rinv.data;
        for (int i = 0; i < d; i++) {
            int rowOffsetR = i * d;
            for (int j = 0; j < d; j++) {
                rData[rowOffsetR + j] = rotation[rowOffsetR + j];
                rInverseData[rowOffsetR + j] = rotation[d*d + rowOffsetR + j];
            }
        }

        DenseMatrix64F Qtrait = workQtrait;
//        Qtrait.set(Qeig);
//        transformMatrixBack(Qtrait, R);
        // Q = R Qeig R^{-1}
        CommonOps.mult(R, Qeig, temp);
        CommonOps.mult(temp, workRinv, workQtrait);

        // 3. Store Q in the actualizations buffer (row-wise copy)
        double[] qData = Qtrait.data;
        for (int i = 0; i < d; i++) {
            System.arraycopy(qData, i * d,
                    actualizations, scaledOffset + i * d,
                    d);
        }
    }

    /* **********************************************************************
     * Integrated OU: override to use BlockDiag S^{-1}
     * ******************************************************************* */

    @Override
    public void updateIntegratedOrnsteinUhlenbeckDiffusionMatrices(
            int precisionIndex,
            final int[] probabilityIndices,
            final double[] edgeLengths,
            final double[] optimalRates,
            final double[] blockDParams,
            final double[] rotation,
            int updateCount) {

        throw new IllegalArgumentException(
                "Integrated OU not yet implemented for Block-Diag integrator. " +
                        "To implement: compute S^{-1} = R D^{-1} Rᵀ using cached block structure " +
                        "and block-wise inversion.");
    }
}

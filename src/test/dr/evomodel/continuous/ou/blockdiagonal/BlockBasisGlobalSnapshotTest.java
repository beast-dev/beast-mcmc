package test.dr.evomodel.continuous.ou.blockdiagonal;

import dr.evomodel.continuous.ou.blockdiagonal.BlockBasisGlobalSnapshot;
import dr.evomodel.continuous.ou.blockdiagonal.BlockDiagonalSelectionMatrixParameterization;
import dr.inference.model.AbstractBlockDiagonalTwoByTwoMatrixParameter;
import dr.inference.model.BlockDiagonalPolarStableMatrixParameter;
import dr.inference.model.MatrixParameter;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Parameter;
import junit.framework.TestCase;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

import java.util.Random;

/**
 * Rotated/block-basis OU work, Phase 1 (see
 * private-code/rotatedOU/rotated_ou_block_basis_plan.md).
 *
 * <p>Checks that {@link BlockBasisGlobalSnapshot}'s D-basis stationary
 * covariance {@code Vtilde = T V T^T} (T = R^{-1}) agrees with two
 * independent computations: (1) a brute-force dense Lyapunov solve plus
 * a dense T V T^T sandwich, using no block-diagonal machinery at all, and
 * (2) the value {@link BlockDiagonalSelectionMatrixParameterization} already
 * computes internally for the identical {@code (R, B, Q)}.</p>
 */
public class BlockBasisGlobalSnapshotTest extends TestCase {

    private static final double TOL = 1e-8;
    private static final int TRIALS = 8;
    private static final int DIM = 3;

    public void testStationaryCovarianceDBasisMatchesBruteForceAndProduction() {
        final Random random = new Random(42);
        for (int trial = 0; trial < TRIALS; ++trial) {
            final MatrixParameter rotationParam = randomWellConditionedRotation(random);
            final Parameter scalarParam = new Parameter.Default(new double[]{-0.5 - random.nextDouble()});
            final Parameter rhoParam = new Parameter.Default(new double[]{-0.3 - 0.6 * random.nextDouble()});
            final Parameter thetaParam = new Parameter.Default(new double[]{(random.nextDouble() - 0.5) * 1.2});
            final Parameter tParam = new Parameter.Default(new double[]{(random.nextDouble() - 0.5) * 0.4});

            final BlockDiagonalPolarStableMatrixParameter blockParameter =
                    new BlockDiagonalPolarStableMatrixParameter(
                            "block" + trial, rotationParam, scalarParam, rhoParam, thetaParam, tParam);

            final MatrixParameter diffusion = randomSpdMatrix(random);

            final BlockBasisGlobalSnapshot snapshot = new BlockBasisGlobalSnapshot(blockParameter, diffusion);
            snapshot.ensureCurrent();

            final DenseMatrix64F vTildeFromSnapshot = snapshot.getStationaryCovarianceDBasis();

            // (1) Brute-force: dense S = R B R^{-1}, solve S V + V S^T + Q = 0 densely,
            // then rotate Vtilde = T V T^T using the snapshot's own T (no block-diagonal solver involved).
            final DenseMatrix64F denseS = denseDriftMatrix(blockParameter);
            final DenseMatrix64F denseQ = toDense(diffusion);
            final DenseMatrix64F denseV = solveDenseLyapunov(denseS, denseQ);
            final DenseMatrix64F vTildeBruteForce = new DenseMatrix64F(DIM, DIM);
            sandwich(snapshot.getRotationInverse(), denseV, vTildeBruteForce);

            assertMatricesClose(
                    "trial " + trial + ": brute-force T V T^T vs snapshot Vtilde",
                    vTildeBruteForce, vTildeFromSnapshot);

            // (2) Cross-check against BlockDiagonalSelectionMatrixParameterization's own internal value.
            final BlockDiagonalSelectionMatrixParameterization production =
                    new BlockDiagonalSelectionMatrixParameterization(blockParameter);
            final double[] productionVTildeFlat = new double[DIM * DIM];
            production.fillStationaryCovarianceDBasisFlatForTesting(diffusion, productionVTildeFlat);
            final DenseMatrix64F vTildeProduction = DenseMatrix64F.wrap(DIM, DIM, productionVTildeFlat);

            assertMatricesClose(
                    "trial " + trial + ": production Vtilde vs snapshot Vtilde",
                    vTildeProduction, vTildeFromSnapshot);

            // Sanity check on the Jacobian term: log|det R| + log|det R^{-1}| == 0.
            final double logAbsDetR = snapshot.getLogAbsDetR();
            final double logAbsDetRinv = Math.log(Math.abs(CommonOps.det(snapshot.getRotationInverse())));
            assertEquals("trial " + trial + ": log|det R| + log|det R^-1| should vanish",
                    0.0, logAbsDetR + logAbsDetRinv, 1e-9);
        }
    }

    private static void assertMatricesClose(final String message,
                                            final DenseMatrix64F expected,
                                            final DenseMatrix64F actual) {
        for (int i = 0; i < DIM; ++i) {
            for (int j = 0; j < DIM; ++j) {
                assertEquals(message + " [" + i + "," + j + "]",
                        expected.get(i, j), actual.get(i, j), TOL);
            }
        }
    }

    private static DenseMatrix64F denseDriftMatrix(final AbstractBlockDiagonalTwoByTwoMatrixParameter blockParameter) {
        final double[] r = new double[DIM * DIM];
        final double[] rinv = new double[DIM * DIM];
        blockParameter.fillRAndRinv(r, rinv);
        final double[] decompositionBlockD = blockParameter.getBlockDiagonalDecomposition().getBlockDiagonal();
        final DenseMatrix64F bMatrix = expandCompressedTridiagonal(decompositionBlockD, DIM);
        final DenseMatrix64F rMatrix = DenseMatrix64F.wrap(DIM, DIM, r);
        final DenseMatrix64F rinvMatrix = DenseMatrix64F.wrap(DIM, DIM, rinv);
        final DenseMatrix64F temp = new DenseMatrix64F(DIM, DIM);
        final DenseMatrix64F out = new DenseMatrix64F(DIM, DIM);
        CommonOps.mult(rMatrix, bMatrix, temp);
        CommonOps.mult(temp, rinvMatrix, out);
        return out;
    }

    private static DenseMatrix64F expandCompressedTridiagonal(final double[] compressed, final int dim) {
        // Layout: [diag(0..dim-1) | upper(0..dim-2) | lower(0..dim-2)], matching
        // BlockDiagonalLyapunovSolver's documented compressed tridiagonal format.
        final DenseMatrix64F out = new DenseMatrix64F(dim, dim);
        for (int i = 0; i < dim; ++i) {
            out.set(i, i, compressed[i]);
        }
        final int upperOffset = dim;
        final int lowerOffset = dim + dim - 1;
        for (int i = 0; i < dim - 1; ++i) {
            out.set(i, i + 1, compressed[upperOffset + i]);
            out.set(i + 1, i, compressed[lowerOffset + i]);
        }
        return out;
    }

    private static DenseMatrix64F solveDenseLyapunov(final DenseMatrix64F s, final DenseMatrix64F q) {
        // Solves S V + V S^T = Q (note: +Q, not -Q), matching
        // BlockDiagonalLyapunovSolver's convention of "D*Sigma + Sigma*D^T = V"
        // where its "V" argument is the RHS (called q/qDBasis at call sites),
        // not the stationary-covariance unknown.
        final int n = s.numRows;
        final int n2 = n * n;
        final DenseMatrix64F bigA = new DenseMatrix64F(n2, n2);
        final DenseMatrix64F bigB = new DenseMatrix64F(n2, 1);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                final int row = i * n + j;
                bigB.set(row, 0, q.get(i, j));
                for (int p = 0; p < n; ++p) {
                    // from (S*V)(i,j) = sum_p S(i,p) V(p,j): coefficient on V(p,j), col = p*n+j
                    bigA.add(row, p * n + j, s.get(i, p));
                    // from (V*S^T)(i,j) = sum_p V(i,p) S(j,p): coefficient on V(i,p), col = i*n+p
                    bigA.add(row, i * n + p, s.get(j, p));
                }
            }
        }
        final DenseMatrix64F bigX = new DenseMatrix64F(n2, 1);
        final boolean solved = CommonOps.solve(bigA, bigB, bigX);
        if (!solved) {
            throw new IllegalStateException("Brute-force dense Lyapunov solve failed (singular system).");
        }
        final DenseMatrix64F vOut = new DenseMatrix64F(n, n);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                vOut.set(i, j, bigX.get(i * n + j, 0));
            }
        }
        return vOut;
    }

    private static void sandwich(final DenseMatrix64F t, final DenseMatrix64F v, final DenseMatrix64F out) {
        final DenseMatrix64F temp = new DenseMatrix64F(t.numRows, v.numCols);
        CommonOps.mult(t, v, temp);
        CommonOps.multTransB(temp, t, out);
    }

    private static DenseMatrix64F toDense(final MatrixParameterInterface parameter) {
        final int dim = parameter.getRowDimension();
        final DenseMatrix64F out = new DenseMatrix64F(dim, dim);
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                out.set(i, j, parameter.getParameterValue(i, j));
            }
        }
        return out;
    }

    private static MatrixParameter randomWellConditionedRotation(final Random random) {
        final MatrixParameter rotation = new MatrixParameter("R", DIM, DIM);
        for (int i = 0; i < DIM; ++i) {
            for (int j = 0; j < DIM; ++j) {
                final double offDiagonal = (i == j) ? 0.0 : (random.nextDouble() - 0.5) * 0.4;
                rotation.setParameterValue(i, j, (i == j ? 1.0 : 0.0) + offDiagonal);
            }
        }
        return rotation;
    }

    private static MatrixParameter randomSpdMatrix(final Random random) {
        final double[][] a = new double[DIM][DIM];
        for (int i = 0; i < DIM; ++i) {
            for (int j = 0; j < DIM; ++j) {
                a[i][j] = random.nextDouble() - 0.5;
            }
        }
        final MatrixParameter q = new MatrixParameter("Q", DIM, DIM);
        for (int i = 0; i < DIM; ++i) {
            for (int j = 0; j < DIM; ++j) {
                double sum = 0.0;
                for (int k = 0; k < DIM; ++k) {
                    sum += a[i][k] * a[j][k];
                }
                if (i == j) {
                    sum += 0.5;
                }
                q.setParameterValue(i, j, sum);
            }
        }
        return q;
    }
}

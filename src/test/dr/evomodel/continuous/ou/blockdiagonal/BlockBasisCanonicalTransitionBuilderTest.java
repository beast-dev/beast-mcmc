package test.dr.evomodel.continuous.ou.blockdiagonal;

import dr.evomodel.continuous.ou.blockdiagonal.BlockBasisCanonicalTransitionBuilder;
import dr.evomodel.continuous.ou.blockdiagonal.BlockBasisGlobalSnapshot;
import dr.evomodel.continuous.ou.blockdiagonal.BlockDiagonalSelectionMatrixParameterization;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.math.GaussianFormConverter;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.message.CanonicalGaussianMessageOps;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.message.CanonicalGaussianState;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.message.CanonicalGaussianTransition;
import dr.inference.model.BlockDiagonalPolarStableMatrixParameter;
import dr.inference.model.MatrixParameter;
import dr.inference.model.Parameter;
import junit.framework.TestCase;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

import java.util.Random;

/**
 * Rotated/block-basis OU work, Phase 2 (see
 * private-code/rotatedOU/rotated_ou_block_basis_plan.md).
 *
 * <p>Builds the canonical OU log-likelihood for a small fixed 3-tip tree
 * {@code ((A:1,B:1):1,C:2)}, fully-observed tips only, two independent ways:
 *
 * <ul>
 *   <li><b>Reference</b>: the existing, unmodified
 *       {@link BlockDiagonalSelectionMatrixParameterization#fillCanonicalTransition}
 *       (natural coordinates, per-branch dense rotation), with the standard
 *       tip-conditioning / message-combine / root-prior-integration primitives
 *       from {@link CanonicalGaussianMessageOps}.</li>
 *   <li><b>New</b>: {@link BlockBasisCanonicalTransitionBuilder} (Ytilde
 *       coordinates, no per-branch rotation), using the same combine/push
 *       primitives, with tip values pre-transformed by {@code T(y - beta)}
 *       and a root prior expressed directly in the D-basis via the cached
 *       {@link BlockBasisGlobalSnapshot}.</li>
 * </ul>
 *
 * <p>Both sides hand-roll the same small postorder message-passing recursion
 * (rather than going through the full tree-likelihood/XML machinery, which
 * Phase 2 deliberately does not touch yet) so the comparison isolates exactly
 * the one thing that changed: how the per-branch transition is built.</p>
 */
public class BlockBasisCanonicalTransitionBuilderTest extends TestCase {

    private static final int DIM = 3;
    private static final double TOL = 1e-7;
    private static final int TRIALS = 5;
    private static final double KAPPA0 = 10.0;

    public void testRotatedBasisLikelihoodMatchesNaturalBasis() {
        final Random random = new Random(7);
        for (int trial = 0; trial < TRIALS; ++trial) {
            final MatrixParameter rotationParam = randomWellConditionedRotation(random);
            final Parameter scalarParam = new Parameter.Default(new double[]{0.5 + random.nextDouble()});
            final Parameter rhoParam = new Parameter.Default(new double[]{0.3 + 0.6 * random.nextDouble()});
            final Parameter thetaParam = new Parameter.Default(new double[]{(random.nextDouble() - 0.5) * 1.2});
            final Parameter tParam = new Parameter.Default(new double[]{(random.nextDouble() - 0.5) * 0.4});
            final BlockDiagonalPolarStableMatrixParameter blockParameter =
                    new BlockDiagonalPolarStableMatrixParameter(
                            "block" + trial, rotationParam, scalarParam, rhoParam, thetaParam, tParam);

            final MatrixParameter diffusion = randomSpdMatrix(random);
            final double[] beta = {random.nextGaussian(), random.nextGaussian(), random.nextGaussian()};
            final double[][] tipValues = {
                    {11.0, 12.0, -2.0},
                    {1.0, 2.0, 5.0},
                    {3.0, 3.0, -4.0}
            };

            final double referenceLogL = computeReferenceLogLikelihood(blockParameter, diffusion, beta, tipValues);
            final double rotatedLogL = computeRotatedLogLikelihood(
                    blockParameter, diffusion, beta, tipValues);

            // Tip values are fixed numbers in natural (Y) coordinates; evaluating their density via the
            // Ytilde-basis machinery requires one log|det T| change-of-variables term per observed tip
            // (T = R^{-1}). This is exact and tip-count-only: the per-branch combine/push Jacobians do not
            // independently survive into the final tree likelihood, they telescope down to this single term.
            final BlockBasisGlobalSnapshot snapshot = new BlockBasisGlobalSnapshot(blockParameter, diffusion);
            snapshot.ensureCurrent();
            final double logAbsDetT = -snapshot.getLogAbsDetR();
            final double expectedRotatedLogL = referenceLogL - tipValues.length * logAbsDetT;

            assertEquals("trial " + trial + ": rotated-basis logL must match natural-basis logL up to the tip-count Jacobian",
                    expectedRotatedLogL, rotatedLogL, TOL);
            assertTrue("trial " + trial + ": logL must be finite", Double.isFinite(referenceLogL));
        }
    }

    /** Natural-basis reference, using the existing unmodified production transition builder. */
    private double computeReferenceLogLikelihood(final BlockDiagonalPolarStableMatrixParameter blockParameter,
                                                  final MatrixParameter diffusion,
                                                  final double[] beta,
                                                  final double[][] tipValues) {
        final BlockDiagonalSelectionMatrixParameterization production =
                new BlockDiagonalSelectionMatrixParameterization(blockParameter);

        final CanonicalGaussianTransition transitionLen1 = new CanonicalGaussianTransition(DIM);
        production.fillCanonicalTransition(diffusion, beta, 1.0, transitionLen1);
        final CanonicalGaussianTransition transitionLen2 = new CanonicalGaussianTransition(DIM);
        production.fillCanonicalTransition(diffusion, beta, 2.0, transitionLen2);

        final CanonicalGaussianMessageOps.Workspace workspace = new CanonicalGaussianMessageOps.Workspace(DIM);

        final CanonicalGaussianState fromA = new CanonicalGaussianState(DIM);
        CanonicalGaussianMessageOps.conditionOnObservedSecondBlock(transitionLen1, tipValues[0], fromA);
        final CanonicalGaussianState fromB = new CanonicalGaussianState(DIM);
        CanonicalGaussianMessageOps.conditionOnObservedSecondBlock(transitionLen1, tipValues[1], fromB);

        final CanonicalGaussianState n1 = new CanonicalGaussianState(DIM);
        CanonicalGaussianMessageOps.combineStates(fromA, fromB, n1);

        final CanonicalGaussianState fromN1AtRoot = new CanonicalGaussianState(DIM);
        CanonicalGaussianMessageOps.pushBackward(n1, transitionLen1, workspace, fromN1AtRoot);

        final CanonicalGaussianState fromC = new CanonicalGaussianState(DIM);
        CanonicalGaussianMessageOps.conditionOnObservedSecondBlock(transitionLen2, tipValues[2], fromC);

        final CanonicalGaussianState rootMessage = new CanonicalGaussianState(DIM);
        CanonicalGaussianMessageOps.combineStates(fromN1AtRoot, fromC, rootMessage);

        final BlockBasisGlobalSnapshot snapshot = new BlockBasisGlobalSnapshot(blockParameter, diffusion);
        snapshot.ensureCurrent();
        final DenseMatrix64F vNatural = naturalStationaryCovariance(snapshot);

        final double[] priorCovFlat = new double[DIM * DIM];
        for (int i = 0; i < DIM * DIM; ++i) {
            priorCovFlat[i] = vNatural.data[i] / KAPPA0;
        }
        final GaussianFormConverter.Workspace converterWorkspace = new GaussianFormConverter.Workspace();
        converterWorkspace.ensureDim(DIM);
        final CanonicalGaussianState priorState = new CanonicalGaussianState(DIM);
        GaussianFormConverter.fillStateFromMoments(beta, priorCovFlat, DIM, converterWorkspace, priorState);

        final CanonicalGaussianState combined = new CanonicalGaussianState(DIM);
        CanonicalGaussianMessageOps.combineStates(rootMessage, priorState, combined);
        return CanonicalGaussianMessageOps.normalizationShift(combined, workspace);
    }

    /** Rotated (Ytilde) basis, using the new Phase 2 block-basis transition builder. */
    private double computeRotatedLogLikelihood(final BlockDiagonalPolarStableMatrixParameter blockParameter,
                                                final MatrixParameter diffusion,
                                                final double[] beta,
                                                final double[][] tipValues) {
        final BlockBasisGlobalSnapshot snapshot = new BlockBasisGlobalSnapshot(blockParameter, diffusion);
        snapshot.ensureCurrent();
        final BlockBasisCanonicalTransitionBuilder builder = new BlockBasisCanonicalTransitionBuilder(
                DIM, snapshot.getBlockStarts(), snapshot.getBlockSizes());

        final CanonicalGaussianTransition transitionLen1 = new CanonicalGaussianTransition(DIM);
        builder.fillCanonicalTransition(snapshot, 1.0, transitionLen1);
        final CanonicalGaussianTransition transitionLen2 = new CanonicalGaussianTransition(DIM);
        builder.fillCanonicalTransition(snapshot, 2.0, transitionLen2);

        final double[][] tipValuesTilde = new double[tipValues.length][DIM];
        for (int t = 0; t < tipValues.length; ++t) {
            transformToTilde(snapshot, tipValues[t], beta, tipValuesTilde[t]);
        }

        final CanonicalGaussianMessageOps.Workspace workspace = new CanonicalGaussianMessageOps.Workspace(DIM);

        final CanonicalGaussianState fromA = new CanonicalGaussianState(DIM);
        CanonicalGaussianMessageOps.conditionOnObservedSecondBlock(transitionLen1, tipValuesTilde[0], fromA);
        final CanonicalGaussianState fromB = new CanonicalGaussianState(DIM);
        CanonicalGaussianMessageOps.conditionOnObservedSecondBlock(transitionLen1, tipValuesTilde[1], fromB);

        final CanonicalGaussianState n1 = new CanonicalGaussianState(DIM);
        CanonicalGaussianMessageOps.combineStates(fromA, fromB, n1);

        final CanonicalGaussianState fromN1AtRoot = new CanonicalGaussianState(DIM);
        CanonicalGaussianMessageOps.pushBackward(n1, transitionLen1, workspace, fromN1AtRoot);

        final CanonicalGaussianState fromC = new CanonicalGaussianState(DIM);
        CanonicalGaussianMessageOps.conditionOnObservedSecondBlock(transitionLen2, tipValuesTilde[2], fromC);

        final CanonicalGaussianState rootMessage = new CanonicalGaussianState(DIM);
        CanonicalGaussianMessageOps.combineStates(fromN1AtRoot, fromC, rootMessage);

        final double[] zeroMean = new double[DIM];
        final double[] priorCovFlat = new double[DIM * DIM];
        final DenseMatrix64F vTilde = snapshot.getStationaryCovarianceDBasis();
        for (int i = 0; i < DIM * DIM; ++i) {
            priorCovFlat[i] = vTilde.data[i] / KAPPA0;
        }
        final GaussianFormConverter.Workspace converterWorkspace = new GaussianFormConverter.Workspace();
        converterWorkspace.ensureDim(DIM);
        final CanonicalGaussianState priorState = new CanonicalGaussianState(DIM);
        GaussianFormConverter.fillStateFromMoments(zeroMean, priorCovFlat, DIM, converterWorkspace, priorState);

        final CanonicalGaussianState combined = new CanonicalGaussianState(DIM);
        CanonicalGaussianMessageOps.combineStates(rootMessage, priorState, combined);
        return CanonicalGaussianMessageOps.normalizationShift(combined, workspace);
    }

    private static void transformToTilde(final BlockBasisGlobalSnapshot snapshot,
                                         final double[] y,
                                         final double[] beta,
                                         final double[] out) {
        final DenseMatrix64F t = snapshot.getRotationInverse();
        for (int i = 0; i < DIM; ++i) {
            double sum = 0.0;
            for (int j = 0; j < DIM; ++j) {
                sum += t.get(i, j) * (y[j] - beta[j]);
            }
            out[i] = sum;
        }
    }

    /** {@code Vnatural = R * Vtilde * R^T}, a one-time (non-hot-path) cross-basis transform for the test's own root prior. */
    private static DenseMatrix64F naturalStationaryCovariance(final BlockBasisGlobalSnapshot snapshot) {
        final DenseMatrix64F r = snapshot.getRotation();
        final DenseMatrix64F vTilde = snapshot.getStationaryCovarianceDBasis();
        final DenseMatrix64F temp = new DenseMatrix64F(DIM, DIM);
        final DenseMatrix64F out = new DenseMatrix64F(DIM, DIM);
        CommonOps.mult(r, vTilde, temp);
        CommonOps.multTransB(temp, r, out);
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

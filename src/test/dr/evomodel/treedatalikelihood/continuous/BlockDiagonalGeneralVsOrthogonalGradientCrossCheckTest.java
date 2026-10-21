package test.dr.evomodel.treedatalikelihood.continuous;

import dr.evomodel.continuous.ou.blockdiagonal.BlockDiagonalSelectionMatrixParameterization;
import dr.evomodel.continuous.ou.blockdiagonal.OrthogonalBlockDiagonalSelectionMatrixParameterization;
import dr.evomodel.continuous.ou.canonical.CanonicalBranchWorkspace;
import dr.evomodel.continuous.ou.canonical.CanonicalPreparedBranchHandle;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.message.CanonicalLocalTransitionAdjoints;
import dr.inference.model.GivensRotationMatrixParameter;
import dr.inference.model.MatrixParameter;
import dr.inference.model.OrthogonalBlockDiagonalPolarStableMatrixParameter;
import dr.inference.model.Parameter;
import junit.framework.TestCase;

/**
 * Cross-checks the general (invertible-R) block-diagonal backend against the
 * orthogonal-specialized backend when both wrap the <em>same</em> orthogonal
 * {@code AbstractBlockDiagonalTwoByTwoMatrixParameter} instance.
 *
 * <p>{@link BlockDiagonalSelectionMatrixParameterization} and
 * {@link OrthogonalBlockDiagonalSelectionMatrixParameterization} compute their
 * transition/covariance sandwich adjoints independently (the latter via
 * {@code OrthogonalBlockCovarianceAdjoint}, using a differently-shaped formula
 * than the former's {@code accumulateSpecializedCovariancePullback}). Both
 * {@code accumulateNativeGradientFrom...Flat} entry points return a dense
 * ({@code dimension x dimension}) adjoint for {@code R}'s entries before the
 * angle-space pullback. That dense adjoint is not unique on the orthogonal
 * manifold, so these tests compare it after the common Givens-angle pullback.
 */
public class BlockDiagonalGeneralVsOrthogonalGradientCrossCheckTest extends TestCase {

    public BlockDiagonalGeneralVsOrthogonalGradientCrossCheckTest(final String name) {
        super(name);
    }

    public void testTransitionDenseRotationGradientMatchesGeneralBackend() {
        final Fixture fixture = buildFixture();
        final OrthogonalBlockDiagonalSelectionMatrixParameterization orthogonalParam =
                new OrthogonalBlockDiagonalSelectionMatrixParameterization(
                        fixture.blockSelection, fixture.rotation);
        final BlockDiagonalSelectionMatrixParameterization generalParam =
                new BlockDiagonalSelectionMatrixParameterization(fixture.blockSelection);

        final double[] dLogL_dF = new double[]{
                 0.10, -0.25,  0.08,  0.03,
                -0.07,  0.18, -0.11,  0.04,
                 0.06, -0.02,  0.21, -0.09,
                -0.04,  0.05, -0.13,  0.17
        };
        final double[] dLogL_df = new double[]{0.14, -0.09, 0.07, -0.03};
        final double[] stationaryMean = new double[]{0.3, -0.4, 0.8, -0.2};
        final double dt = 0.37;

        final double[] orthogonalCompressed = new double[orthogonalParam.getCompressedSelectionGradientDimension()];
        final double[] orthogonalRotation = new double[16];
        orthogonalParam.accumulateNativeGradientFromTransitionFlat(
                dt, stationaryMean, dLogL_dF, dLogL_df, orthogonalCompressed, orthogonalRotation);

        final double[] generalCompressed = new double[generalParam.getCompressedSelectionGradientDimension()];
        final double[] generalRotation = new double[16];
        generalParam.accumulateNativeGradientFromTransitionFlat(
                dt, stationaryMean, dLogL_dF, dLogL_df, generalCompressed, generalRotation);

        assertVectorEquals("compressed transition gradient", generalCompressed, orthogonalCompressed, 1.0e-10);

        final double[] generalAngles = new double[fixture.angles.getDimension()];
        fixture.rotation.fillPullBackGradientFlat(generalRotation, 4, generalAngles);
        final double[] orthogonalAngles = new double[fixture.angles.getDimension()];
        fixture.rotation.fillPullBackGradientFlat(orthogonalRotation, 4, orthogonalAngles);
        assertVectorEquals("angle-space transition gradient", generalAngles, orthogonalAngles, 1.0e-10);
    }

    public void testCovarianceDenseRotationGradientMatchesGeneralBackend() {
        final Fixture fixture = buildFixture();
        final OrthogonalBlockDiagonalSelectionMatrixParameterization orthogonalParam =
                new OrthogonalBlockDiagonalSelectionMatrixParameterization(
                        fixture.blockSelection, fixture.rotation);
        final BlockDiagonalSelectionMatrixParameterization generalParam =
                new BlockDiagonalSelectionMatrixParameterization(fixture.blockSelection);

        final MatrixParameter diffusionMatrix = buildDiffusionMatrix();
        final double[] dLogL_dV = new double[]{
                0.40, -0.06,  0.03,  0.01,
               -0.06,  0.28, -0.04,  0.02,
                0.03, -0.04,  0.35, -0.07,
                0.01,  0.02, -0.07,  0.31
        };
        final double dt = 0.29;

        final double[] orthogonalCompressed = new double[orthogonalParam.getCompressedSelectionGradientDimension()];
        final double[] orthogonalRotation = new double[16];
        orthogonalParam.accumulateNativeGradientFromCovarianceStationaryFlat(
                diffusionMatrix, dt, dLogL_dV, orthogonalCompressed, orthogonalRotation);

        final double[] generalCompressed = new double[generalParam.getCompressedSelectionGradientDimension()];
        final double[] generalRotation = new double[16];
        generalParam.accumulateNativeGradientFromCovarianceStationaryFlat(
                diffusionMatrix, dt, dLogL_dV, generalCompressed, generalRotation);

        assertVectorEquals("compressed covariance gradient", generalCompressed, orthogonalCompressed, 1.0e-8);

        final double[] generalAngles = new double[fixture.angles.getDimension()];
        fixture.rotation.fillPullBackGradientFlat(generalRotation, 4, generalAngles);
        final double[] orthogonalAngles = new double[fixture.angles.getDimension()];
        fixture.rotation.fillPullBackGradientFlat(orthogonalRotation, 4, orthogonalAngles);
        assertVectorEquals("angle-space covariance gradient", generalAngles, orthogonalAngles, 1.0e-10);
    }

    private static Fixture buildFixture() {
        final Parameter angles = new Parameter.Default("angles.crossCheck",
                new double[]{0.10, -0.12, 0.08, 0.04, -0.06, 0.11});
        final GivensRotationMatrixParameter rotation =
                new GivensRotationMatrixParameter("rotation.crossCheck", 4, angles);
        final Parameter scalar = new Parameter.Default("scalar.crossCheck", new double[0]);
        final Parameter rho = new Parameter.Default("rho.crossCheck", new double[]{0.85, 1.10});
        final Parameter theta = new Parameter.Default("theta.crossCheck", new double[]{0.55, -0.35});
        final Parameter t = new Parameter.Default("t.crossCheck", new double[]{0.27, -0.19});

        final OrthogonalBlockDiagonalPolarStableMatrixParameter blockSelection =
                new OrthogonalBlockDiagonalPolarStableMatrixParameter(
                        "selection.crossCheck", rotation, scalar, rho, theta, t);

        return new Fixture(blockSelection, rotation, angles);
    }

    private static final class Fixture {
        private final OrthogonalBlockDiagonalPolarStableMatrixParameter blockSelection;
        private final GivensRotationMatrixParameter rotation;
        private final Parameter angles;

        private Fixture(final OrthogonalBlockDiagonalPolarStableMatrixParameter blockSelection,
                        final GivensRotationMatrixParameter rotation,
                        final Parameter angles) {
            this.blockSelection = blockSelection;
            this.rotation = rotation;
            this.angles = angles;
        }
    }

    /**
     * Adjudicates the covariance-pullback discrepancy found by
     * {@link #testCovarianceDenseRotationGradientMatchesGeneralBackend} against an
     * independent finite-difference ground truth. Defines a synthetic scalar loss
     * {@code L(theta) = dLogL_dV . V(R(theta))} (V = transition covariance, shared/
     * unambiguous machinery for both backends), computes dL/dtheta by perturbing the
     * Givens angles directly, and compares against each backend's dense R-gradient
     * pulled back through the same {@code OrthogonalMatrixProvider.fillPullBackGradientFlat}
     * Jacobian.
     */
    public void testCovarianceRotationGradientAgainstFiniteDifference() {
        final Fixture fixture = buildFixture();
        final MatrixParameter diffusionMatrix = buildDiffusionMatrix();
        final double dt = 0.29;
        final double[] dLogL_dV = new double[]{
                0.40, -0.06,  0.03,  0.01,
               -0.06,  0.28, -0.04,  0.02,
                0.03, -0.04,  0.35, -0.07,
                0.01,  0.02, -0.07,  0.31
        };
        final int dimension = 4;
        final int angleCount = fixture.angles.getDimension();

        final BlockDiagonalSelectionMatrixParameterization refParam =
                new BlockDiagonalSelectionMatrixParameterization(fixture.blockSelection);
        final double[] originalAngles = new double[angleCount];
        for (int i = 0; i < angleCount; ++i) {
            originalAngles[i] = fixture.angles.getParameterValue(i);
        }

        final double eps = 1.0e-6;
        final double[] fdAngleGradient = new double[angleCount];
        final double[] vFlat = new double[dimension * dimension];
        for (int k = 0; k < angleCount; ++k) {
            fixture.angles.setParameterValue(k, originalAngles[k] + eps);
            refParam.fillTransitionCovarianceFlat(diffusionMatrix, dt, vFlat);
            final double lPlus = dot(dLogL_dV, vFlat);

            fixture.angles.setParameterValue(k, originalAngles[k] - eps);
            refParam.fillTransitionCovarianceFlat(diffusionMatrix, dt, vFlat);
            final double lMinus = dot(dLogL_dV, vFlat);

            fixture.angles.setParameterValue(k, originalAngles[k]);
            fdAngleGradient[k] = (lPlus - lMinus) / (2.0 * eps);
        }

        final OrthogonalBlockDiagonalSelectionMatrixParameterization orthogonalParam =
                new OrthogonalBlockDiagonalSelectionMatrixParameterization(
                        fixture.blockSelection, fixture.rotation);
        final BlockDiagonalSelectionMatrixParameterization generalParam =
                new BlockDiagonalSelectionMatrixParameterization(fixture.blockSelection);

        final double[] orthogonalCompressed = new double[orthogonalParam.getCompressedSelectionGradientDimension()];
        final double[] orthogonalRotationDense = new double[dimension * dimension];
        orthogonalParam.accumulateNativeGradientFromCovarianceStationaryFlat(
                diffusionMatrix, dt, dLogL_dV, orthogonalCompressed, orthogonalRotationDense);

        final double[] generalCompressed = new double[generalParam.getCompressedSelectionGradientDimension()];
        final double[] generalRotationDense = new double[dimension * dimension];
        generalParam.accumulateNativeGradientFromCovarianceStationaryFlat(
                diffusionMatrix, dt, dLogL_dV, generalCompressed, generalRotationDense);

        final double[] orthogonalAngleGradient = new double[angleCount];
        fixture.rotation.fillPullBackGradientFlat(orthogonalRotationDense, dimension, orthogonalAngleGradient);
        final double[] generalAngleGradient = new double[angleCount];
        fixture.rotation.fillPullBackGradientFlat(generalRotationDense, dimension, generalAngleGradient);

        System.out.println("FD angle gradient:         " + java.util.Arrays.toString(fdAngleGradient));
        System.out.println("Orthogonal-formula pullback: " + java.util.Arrays.toString(orthogonalAngleGradient));
        System.out.println("General-formula pullback:    " + java.util.Arrays.toString(generalAngleGradient));
        assertVectorEquals("orthogonal covariance angle-gradient vs FD",
                fdAngleGradient, orthogonalAngleGradient, 1.0e-5);
        assertVectorEquals("general covariance angle-gradient vs FD",
                fdAngleGradient, generalAngleGradient, 1.0e-5);
    }

    /**
     * Verifies the "prepared" flavor of the simple covariance/R-gradient pullback
     * ({@code accumulateNativeGradientFromAdjointsPreparedFlat}) separately from the
     * "current" flavor already checked above. This entry point is called by production
     * code ({@code CanonicalPreparedBranchGradientAccumulator}) but, unlike the
     * diffusion-combined variant, is not exercised by any existing orthogonal test —
     * exactly the kind of gap where the fillNativeTransitionCovarianceMatrix field-mismatch
     * bug hid undetected.
     */
    public void testPreparedCovarianceGradientMatchesGeneralBackendAndFiniteDifference() {
        final Fixture fixture = buildFixture();
        final BlockDiagonalSelectionMatrixParameterization refParam =
                new BlockDiagonalSelectionMatrixParameterization(fixture.blockSelection);
        final OrthogonalBlockDiagonalSelectionMatrixParameterization orthogonalParam =
                new OrthogonalBlockDiagonalSelectionMatrixParameterization(fixture.blockSelection, fixture.rotation);
        final BlockDiagonalSelectionMatrixParameterization generalParam =
                new BlockDiagonalSelectionMatrixParameterization(fixture.blockSelection);

        final MatrixParameter diffusionMatrix = buildDiffusionMatrix();
        final double dt = 0.29;
        final int dimension = 4;
        final double[] stationaryMean = new double[]{0.25, -0.35, 0.55, -0.15};
        final double[] dLogL_dOmega = new double[]{
                0.40, -0.06,  0.03,  0.01,
               -0.06,  0.28, -0.04,  0.02,
                0.03, -0.04,  0.35, -0.07,
                0.01,  0.02, -0.07,  0.31
        };

        final double[] orthogonalCompressed = new double[orthogonalParam.getCompressedSelectionGradientDimension()];
        final double[] orthogonalRotationDense = new double[dimension * dimension];
        accumulatePreparedCovarianceGradient(
                orthogonalParam, diffusionMatrix, stationaryMean, dLogL_dOmega, dt, orthogonalCompressed, orthogonalRotationDense);

        final double[] generalCompressed = new double[generalParam.getCompressedSelectionGradientDimension()];
        final double[] generalRotationDense = new double[dimension * dimension];
        accumulatePreparedCovarianceGradient(
                generalParam, diffusionMatrix, stationaryMean, dLogL_dOmega, dt, generalCompressed, generalRotationDense);

        assertVectorEquals("prepared compressed covariance gradient", generalCompressed, orthogonalCompressed, 1.0e-10);

        final double[] generalAngleFromPrepared = new double[fixture.angles.getDimension()];
        fixture.rotation.fillPullBackGradientFlat(generalRotationDense, dimension, generalAngleFromPrepared);
        final double[] orthogonalAngleFromPrepared = new double[fixture.angles.getDimension()];
        fixture.rotation.fillPullBackGradientFlat(orthogonalRotationDense, dimension, orthogonalAngleFromPrepared);
        assertVectorEquals("prepared angle-space covariance gradient",
                generalAngleFromPrepared, orthogonalAngleFromPrepared, 1.0e-10);

        // Cross-check against finite differences too, not just against each other.
        final int angleCount = fixture.angles.getDimension();
        final double[] originalAngles = new double[angleCount];
        for (int i = 0; i < angleCount; ++i) {
            originalAngles[i] = fixture.angles.getParameterValue(i);
        }
        final double eps = 1.0e-6;
        final double[] fdAngleGradient = new double[angleCount];
        final double[] vFlat = new double[dimension * dimension];
        for (int k = 0; k < angleCount; ++k) {
            fixture.angles.setParameterValue(k, originalAngles[k] + eps);
            refParam.fillTransitionCovarianceFlat(diffusionMatrix, dt, vFlat);
            final double lPlus = dot(dLogL_dOmega, vFlat);

            fixture.angles.setParameterValue(k, originalAngles[k] - eps);
            refParam.fillTransitionCovarianceFlat(diffusionMatrix, dt, vFlat);
            final double lMinus = dot(dLogL_dOmega, vFlat);

            fixture.angles.setParameterValue(k, originalAngles[k]);
            fdAngleGradient[k] = (lPlus - lMinus) / (2.0 * eps);
        }
        final double[] orthogonalAngleGradient = new double[angleCount];
        fixture.rotation.fillPullBackGradientFlat(orthogonalRotationDense, dimension, orthogonalAngleGradient);
        assertVectorEquals("prepared angle-space covariance gradient vs FD", fdAngleGradient, orthogonalAngleGradient, 1.0e-5);
    }

    private static void accumulatePreparedCovarianceGradient(final BlockDiagonalSelectionMatrixParameterization parameterization,
                                                             final MatrixParameter diffusionMatrix,
                                                             final double[] stationaryMean,
                                                             final double[] dLogL_dOmega,
                                                             final double dt,
                                                             final double[] compressedOut,
                                                             final double[] rotationOut) {
        final CanonicalPreparedBranchHandle prepared = parameterization.createPreparedBranchHandle();
        final CanonicalBranchWorkspace workspace = parameterization.createBranchWorkspace();
        parameterization.prepareBranch(dt, stationaryMean, prepared);

        final CanonicalLocalTransitionAdjoints adjoints = new CanonicalLocalTransitionAdjoints(4);
        System.arraycopy(dLogL_dOmega, 0, adjoints.dLogL_dOmega, 0, dLogL_dOmega.length);

        parameterization.accumulateNativeGradientFromAdjointsPreparedFlat(
                prepared, diffusionMatrix, adjoints, workspace, compressedOut, rotationOut);
    }

    /**
     * Verifies {@code accumulatePreparedCovarianceAndDiffusionGradient} — the one
     * remaining, un-merged override that jointly produces the R-gradient AND the
     * diffusion (Q) gradient in one call. The orthogonal formula reuses the yAdjoint
     * computed inside its own covariance-pullback step for the diffusion gradient;
     * the general formula recomputes it via a separate fillSymmetricSandwichOutputAdjoint
     * + fillDiffusionGradientDBasisPrepared call. Cross-checked against finite
     * differences on Q, since (unlike R) Q has no manifold-constraint gauge freedom.
     */
    public void testPreparedDiffusionCombinedGradientAgainstFiniteDifference() {
        final Fixture fixture = buildFixture();
        final BlockDiagonalSelectionMatrixParameterization refParam =
                new BlockDiagonalSelectionMatrixParameterization(fixture.blockSelection);
        final OrthogonalBlockDiagonalSelectionMatrixParameterization orthogonalParam =
                new OrthogonalBlockDiagonalSelectionMatrixParameterization(fixture.blockSelection, fixture.rotation);
        final BlockDiagonalSelectionMatrixParameterization generalParam =
                new BlockDiagonalSelectionMatrixParameterization(fixture.blockSelection);

        final MatrixParameter diffusionMatrix = buildDiffusionMatrix();
        final double dt = 0.29;
        final int dimension = 4;
        final double[] stationaryMean = new double[]{0.25, -0.35, 0.55, -0.15};
        final double[] dLogL_dOmega = new double[]{
                0.40, -0.06,  0.03,  0.01,
               -0.06,  0.28, -0.04,  0.02,
                0.03, -0.04,  0.35, -0.07,
                0.01,  0.02, -0.07,  0.31
        };

        final double[] orthogonalCompressed = new double[orthogonalParam.getCompressedSelectionGradientDimension()];
        final double[] orthogonalRotationDense = new double[dimension * dimension];
        final double[] orthogonalDiffusionGrad = new double[dimension * dimension];
        accumulatePreparedCovarianceAndDiffusionGradient(
                orthogonalParam, diffusionMatrix, stationaryMean, dLogL_dOmega, dt,
                orthogonalCompressed, orthogonalRotationDense, orthogonalDiffusionGrad);

        final double[] generalCompressed = new double[generalParam.getCompressedSelectionGradientDimension()];
        final double[] generalRotationDense = new double[dimension * dimension];
        final double[] generalDiffusionGrad = new double[dimension * dimension];
        accumulatePreparedCovarianceAndDiffusionGradient(
                generalParam, diffusionMatrix, stationaryMean, dLogL_dOmega, dt,
                generalCompressed, generalRotationDense, generalDiffusionGrad);

        assertVectorEquals("prepared+diffusion compressed gradient", generalCompressed, orthogonalCompressed, 1.0e-10);
        // Note: unlike the already-merged methods, the dense rotation gradient is NOT
        // expected to match bit-for-bit here yet (this method still uses the old
        // orthogonal-specific formula) — only the angle-space projection below must agree.

        // R-gradient: cross-check against FD via the angle pullback (gauge-safe).
        final int angleCount = fixture.angles.getDimension();
        final double[] originalAngles = new double[angleCount];
        for (int i = 0; i < angleCount; ++i) {
            originalAngles[i] = fixture.angles.getParameterValue(i);
        }
        final double eps = 1.0e-6;
        final double[] fdAngleGradient = new double[angleCount];
        final double[] vFlat = new double[dimension * dimension];
        for (int k = 0; k < angleCount; ++k) {
            fixture.angles.setParameterValue(k, originalAngles[k] + eps);
            refParam.fillTransitionCovarianceFlat(diffusionMatrix, dt, vFlat);
            final double lPlus = dot(dLogL_dOmega, vFlat);

            fixture.angles.setParameterValue(k, originalAngles[k] - eps);
            refParam.fillTransitionCovarianceFlat(diffusionMatrix, dt, vFlat);
            final double lMinus = dot(dLogL_dOmega, vFlat);

            fixture.angles.setParameterValue(k, originalAngles[k]);
            fdAngleGradient[k] = (lPlus - lMinus) / (2.0 * eps);
        }
        final double[] orthogonalAngleGradient = new double[angleCount];
        fixture.rotation.fillPullBackGradientFlat(orthogonalRotationDense, dimension, orthogonalAngleGradient);
        assertVectorEquals("prepared+diffusion angle-space R-gradient vs FD",
                fdAngleGradient, orthogonalAngleGradient, 1.0e-5);

        // Diffusion (Q) gradient: cross-check via random symmetric directional derivatives
        // (no gauge freedom for Q, so this must match directly, not just after a projection).
        final java.util.Random random = new java.util.Random(7);
        for (int trial = 0; trial < 4; ++trial) {
            final double[] direction = randomSymmetricDirection(dimension, random);
            final double lPlus = evaluateTransitionCovarianceLoss(
                    refParam, diffusionMatrix, dt, dLogL_dOmega, direction, eps);
            final double lMinus = evaluateTransitionCovarianceLoss(
                    refParam, diffusionMatrix, dt, dLogL_dOmega, direction, -eps);
            final double fd = (lPlus - lMinus) / (2.0 * eps);

            final double orthogonalDirectional = dot(orthogonalDiffusionGrad, direction);
            final double generalDirectional = dot(generalDiffusionGrad, direction);
            final double tolerance = 1.0e-5 * Math.max(1.0, Math.abs(fd));
            assertEquals("prepared+diffusion orthogonal Q-gradient directional, trial " + trial,
                    fd, orthogonalDirectional, tolerance);
            assertEquals("prepared+diffusion general Q-gradient directional, trial " + trial,
                    fd, generalDirectional, tolerance);
        }
    }

    private static void accumulatePreparedCovarianceAndDiffusionGradient(
            final BlockDiagonalSelectionMatrixParameterization parameterization,
            final MatrixParameter diffusionMatrix,
            final double[] stationaryMean,
            final double[] dLogL_dOmega,
            final double dt,
            final double[] compressedOut,
            final double[] rotationOut,
            final double[] diffusionOut) {
        final CanonicalPreparedBranchHandle prepared = parameterization.createPreparedBranchHandle();
        final CanonicalBranchWorkspace workspace = parameterization.createBranchWorkspace();
        parameterization.prepareBranch(dt, stationaryMean, prepared);

        final CanonicalLocalTransitionAdjoints adjoints = new CanonicalLocalTransitionAdjoints(4);
        System.arraycopy(dLogL_dOmega, 0, adjoints.dLogL_dOmega, 0, dLogL_dOmega.length);

        parameterization.accumulateNativeSelectionAndDiffusionGradientFromAdjointsPreparedFlat(
                prepared, diffusionMatrix, adjoints, workspace, compressedOut, rotationOut, false, diffusionOut);
    }

    private static double dot(final double[] a, final double[] b) {
        double sum = 0.0;
        for (int i = 0; i < a.length; ++i) {
            sum += a[i] * b[i];
        }
        return sum;
    }

    /**
     * Adjudicates the diffusion (Q) gradient formulas between backends against an
     * independent finite-difference ground truth. Unlike the R-gradient case, Q has
     * no manifold-constraint gauge freedom to hide behind: if the two backends'
     * dense dLogL_dQ formulas disagree at all here, one of them is simply wrong.
     * Uses directional derivatives along random symmetric perturbations (rather than
     * per-entry finite differences) so no assumption about either formula's exact
     * indexing/symmetrization convention is needed — only linearity of the gradient.
     */
    public void testDiffusionGradientAgainstFiniteDifference() {
        final Fixture fixture = buildFixture();
        final BlockDiagonalSelectionMatrixParameterization refParam =
                new BlockDiagonalSelectionMatrixParameterization(fixture.blockSelection);
        final OrthogonalBlockDiagonalSelectionMatrixParameterization orthogonalParam =
                new OrthogonalBlockDiagonalSelectionMatrixParameterization(fixture.blockSelection, fixture.rotation);
        final BlockDiagonalSelectionMatrixParameterization generalParam =
                new BlockDiagonalSelectionMatrixParameterization(fixture.blockSelection);

        final MatrixParameter diffusionMatrix = buildDiffusionMatrix();
        final double dt = 0.29;
        final int dimension = 4;
        final double[] dLogL_dV = new double[]{
                0.40, -0.06,  0.03,  0.01,
               -0.06,  0.28, -0.04,  0.02,
                0.03, -0.04,  0.35, -0.07,
                0.01,  0.02, -0.07,  0.31
        };

        final double[] orthogonalDiffusionGrad = new double[dimension * dimension];
        orthogonalParam.accumulateDiffusionGradientFlat(diffusionMatrix, dt, dLogL_dV, false, orthogonalDiffusionGrad);

        final double[] generalDiffusionGrad = new double[dimension * dimension];
        generalParam.accumulateDiffusionGradientFlat(diffusionMatrix, dt, dLogL_dV, false, generalDiffusionGrad);

        final double eps = 1.0e-6;
        final java.util.Random random = new java.util.Random(42);
        for (int trial = 0; trial < 4; ++trial) {
            final double[] direction = randomSymmetricDirection(dimension, random);

            final double lPlus = evaluateTransitionCovarianceLoss(
                    refParam, diffusionMatrix, dt, dLogL_dV, direction, eps);
            final double lMinus = evaluateTransitionCovarianceLoss(
                    refParam, diffusionMatrix, dt, dLogL_dV, direction, -eps);
            final double fd = (lPlus - lMinus) / (2.0 * eps);

            final double orthogonalDirectional = dot(orthogonalDiffusionGrad, direction);
            final double generalDirectional = dot(generalDiffusionGrad, direction);

            final double tolerance = 1.0e-5 * Math.max(1.0, Math.abs(fd));
            assertEquals("orthogonal diffusion directional derivative, trial " + trial,
                    fd, orthogonalDirectional, tolerance);
            assertEquals("general diffusion directional derivative, trial " + trial,
                    fd, generalDirectional, tolerance);
        }
    }

    private static double[] randomSymmetricDirection(final int dimension, final java.util.Random random) {
        final double[] direction = new double[dimension * dimension];
        for (int i = 0; i < dimension; ++i) {
            for (int j = i; j < dimension; ++j) {
                final double value = random.nextGaussian();
                direction[i * dimension + j] = value;
                direction[j * dimension + i] = value;
            }
        }
        return direction;
    }

    private static double evaluateTransitionCovarianceLoss(final BlockDiagonalSelectionMatrixParameterization refParam,
                                                            final MatrixParameter diffusionMatrix,
                                                            final double dt,
                                                            final double[] dLogL_dV,
                                                            final double[] direction,
                                                            final double eps) {
        final int dimension = 4;
        final double[][] original = new double[dimension][dimension];
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                original[i][j] = diffusionMatrix.getParameterValue(i, j);
                diffusionMatrix.setParameterValueQuietly(i, j, original[i][j] + eps * direction[i * dimension + j]);
            }
        }
        diffusionMatrix.fireParameterChangedEvent();

        final double[] vFlat = new double[dimension * dimension];
        refParam.fillTransitionCovarianceFlat(diffusionMatrix, dt, vFlat);
        final double loss = dot(dLogL_dV, vFlat);

        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                diffusionMatrix.setParameterValueQuietly(i, j, original[i][j]);
            }
        }
        diffusionMatrix.fireParameterChangedEvent();

        return loss;
    }

    private static MatrixParameter buildDiffusionMatrix() {
        final MatrixParameter q = new MatrixParameter("q.crossCheck", 4, 4);
        final double[][] values = new double[][]{
                {2.0, 0.20, 0.05, 0.01},
                {0.20, 1.7, 0.12, -0.03},
                {0.05, 0.12, 1.9, 0.18},
                {0.01, -0.03, 0.18, 1.6}
        };
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                q.setParameterValueQuietly(i, j, values[i][j]);
            }
        }
        q.fireParameterChangedEvent();
        return q;
    }

    private static void assertVectorEquals(final String label,
                                           final double[] expected,
                                           final double[] actual,
                                           final double tolerance) {
        assertEquals(label + " length", expected.length, actual.length);
        for (int i = 0; i < expected.length; ++i) {
            assertEquals(label + "[" + i + "]", expected[i], actual[i], tolerance);
        }
    }
}

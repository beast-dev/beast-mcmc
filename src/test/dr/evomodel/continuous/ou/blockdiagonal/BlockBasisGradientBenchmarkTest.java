package test.dr.evomodel.continuous.ou.blockdiagonal;

import dr.evomodel.continuous.ou.blockdiagonal.BlockBasisCanonicalAdjointBuilder;
import dr.evomodel.continuous.ou.blockdiagonal.BlockBasisCanonicalTransitionBuilder;
import dr.evomodel.continuous.ou.blockdiagonal.BlockBasisGlobalSnapshot;
import dr.evomodel.continuous.ou.blockdiagonal.BlockDiagonalSelectionMatrixParameterization;
import dr.evomodel.continuous.canonical.math.GaussianFormConverter;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.message.CanonicalBranchMessageContribution;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.message.CanonicalBranchMessageContributionUtils;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.message.CanonicalGaussianMessageOps;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.message.CanonicalGaussianState;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.message.CanonicalGaussianTransition;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.message.CanonicalLocalTransitionAdjoints;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.message.CanonicalTransitionAdjointUtils;
import dr.inference.model.BlockDiagonalPolarStableMatrixParameter;
import dr.inference.model.MatrixParameter;
import dr.inference.model.Parameter;
import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.Random;

/**
 * Rotated/block-basis OU work, post-Phase-3 benchmark (see
 * private-code/rotatedOU/rotated_ou_block_basis_plan.md).
 *
 * <p>Standalone (not part of {@code ant junit-canonical}) wall-clock benchmark of one
 * <b>full gradient evaluation</b> under the new block-basis backend — forward
 * (post-order) message passing, backward (pre-order/smoothing) message passing, and
 * every pullback (per-branch Frechet/Lyapunov adjoint, degenerate-tip Fisher
 * contribution + tip-transform R-gradient, root-prior Vtilde-gradient, the global
 * finalize/Lyapunov-adjoint solve, and the explicit Jacobian-term gradient) — at trait
 * dimension {@code p = 4, 8, 16} on a balanced 50-tip synthetic tree.
 *
 * <p>This hand-rolls the same post-order/pre-order recursion that
 * {@code BlockBasisCanonicalAdjointBuilderTest} uses for its fixed 3-tip tree (see that
 * file and {@code BlockBasisCanonicalTransitionBuilderTest} for the derivation), just
 * generalized to an arbitrary balanced binary tree, since Phase 4/5 (generalized tip
 * observation model, full tree-likelihood/XML wiring) have not happened yet.
 *
 * <p>Each measured iteration perturbs the rotation and native block parameters by a
 * tiny, monotonically-changing amount before evaluating, so {@link
 * BlockBasisGlobalSnapshot#ensureCurrent()} is forced to redo its work every time
 * (matching the realistic cost of one gradient evaluation at a freshly-proposed
 * parameter point during MCMC/HMC, rather than measuring the cheap no-op it takes when
 * called repeatedly at an unchanged point).
 */
public final class BlockBasisGradientBenchmarkTest extends TestCase {

    private static final int TIP_COUNT = 100;
    private static final int[] DIMENSIONS = {4, 8, 16};
    private static final int WARMUP_ITERS = 200;
    private static final int MEASURE_ITERS = 1000;
    private static final double KAPPA0 = 10.0;
    private static final long SEED = 42L;
    private static final double FD_DELTA = 1.0e-6;
    private static final double RELATIVE_TOL = 5.0e-4;

    public BlockBasisGradientBenchmarkTest(final String name) {
        super(name);
    }

    public static Test suite() {
        return new TestSuite(BlockBasisGradientBenchmarkTest.class);
    }

    public static void main(final String[] args) {
        junit.textui.TestRunner.run(suite());
    }

    public void testBenchmarkReportIsFiniteAndConsistent() {
        final StringBuilder report = new StringBuilder();
        for (final int dim : DIMENSIONS) {
            report.append(runBenchmark(dim));
        }
        final String text = report.toString();
        System.out.print(text);
        assertFalse(text.contains("NaN"));
        assertFalse(text.contains("Infinity"));
        for (final int dim : DIMENSIONS) {
            assertTrue(text.contains("p=" + dim));
        }
    }

    private String runBenchmark(final int dim) {
        final Random random = new Random(SEED + dim);
        final Fixture fixture = buildFixture(dim, random);
        final Node root = buildBalancedTree(random);
        allocateNodeStates(root, dim, true);
        final Scratch scratch = new Scratch(dim, fixture);

        validateCorrectness(root, fixture, scratch);

        double logL = Double.NaN;
        for (int i = 0; i < WARMUP_ITERS; ++i) {
            perturbParameters(i, fixture);
            logL = computeFullGradient(root, fixture, scratch);
        }
        if (!Double.isFinite(logL)) {
            throw new AssertionError("Non-finite benchmark log-likelihood at p=" + dim + ": " + logL);
        }

        final long start = System.nanoTime();
        for (int i = 0; i < MEASURE_ITERS; ++i) {
            perturbParameters(WARMUP_ITERS + i, fixture);
            logL = computeFullGradient(root, fixture, scratch);
        }
        final long elapsedNs = System.nanoTime() - start;
        final double nsPerGradient = elapsedNs / (double) MEASURE_ITERS;

        if (!isFinite(scratch.rGradAcc.data) || !isFinite(scratch.qGradFlatAcc) || !isFinite(scratch.nativeGradAcc)) {
            throw new AssertionError("Non-finite benchmark gradient at p=" + dim);
        }

        final int edgeCount = 2 * TIP_COUNT - 2;
        final StringBuilder sb = new StringBuilder();
        sb.append("Case p=").append(dim).append(", ").append(TIP_COUNT).append(" tips, ")
                .append(edgeCount).append(" edges\n");
        sb.append(String.format("  logL: %.6f%n", logL));
        sb.append(String.format("  full gradient (post-order + pre-order + all pullbacks): %.1f ns (%.3f us)%n",
                nsPerGradient, nsPerGradient / 1000.0));
        sb.append(String.format("  per edge: %.1f ns%n", nsPerGradient / edgeCount));
        sb.append(String.format("  throughput: %.1f gradients/sec%n", 1.0e9 / nsPerGradient));
        sb.append('\n');
        return sb.toString();
    }

    // ------------------------------------------------------------------
    // One full gradient evaluation: forward, backward, all pullbacks.
    // ------------------------------------------------------------------

    private double computeFullGradient(final Node root, final Fixture fixture, final Scratch s) {
        s.snapshot.ensureCurrent();

        for (int t = 0; t < TIP_COUNT; ++t) {
            transformToTilde(s.snapshot, fixture.tipValuesNatural[t], fixture.beta, s.dim, s.tipValuesTilde[t]);
        }

        Arrays.fill(s.compressedDAcc, 0.0);
        CommonOps.fill(s.gSAcc, 0.0);
        CommonOps.fill(s.rGradAcc, 0.0);
        Arrays.fill(s.qGradFlatAcc, 0.0);

        computePostOrder(root, s);

        final DenseMatrix64F vTilde = s.snapshot.getStationaryCovarianceDBasis();
        for (int i = 0; i < s.dim * s.dim; ++i) {
            s.priorCovFlat[i] = vTilde.data[i] / KAPPA0;
        }
        GaussianFormConverter.fillStateFromMoments(s.zeroMean, s.priorCovFlat, s.dim, s.gfcWs, s.priorState);

        CanonicalGaussianMessageOps.combineStates(root.selfMessage, s.priorState, root.nodeFullPosterior);
        final double logL = CanonicalGaussianMessageOps.normalizationShift(root.nodeFullPosterior, s.ws);

        computeBackwardAndPullback(root, s.priorState, s);
        accumulateRootPriorVtildePullback(root.nodeFullPosterior, s.priorState, s);

        s.adjointBuilder.finalizeGlobalPullback(s.snapshot, s.gSAcc, s.compressedDAcc, s.rGradAcc, s.qGradFlatAcc);
        s.adjointBuilder.addJacobianTermGradient(s.snapshot, TIP_COUNT, s.rGradAcc);

        s.blockParameter.chainGradient(s.compressedDAcc, s.nativeGradAcc);
        return logL;
    }

    // ------------------------------------------------------------------
    // Correctness validation, at the fixture's base (unperturbed) parameters, before
    // trusting the timing numbers: the generalized N-tip recursion above is new (Phase
    // 2/3's hand-rolled tests only cover a fixed 3-tip tree), so this checks (a) the
    // rotated-basis logL against an independently-recursed dense/natural-basis
    // reference (validates the forward post-order/pre-order tree traversal generalizes
    // correctly), and (b) a few finite-difference spot checks of the analytic gradient
    // (rho of the first 2x2 block, one R entry, one Q entry) against that same
    // reference (validates the backward pass + pullback machinery).
    // ------------------------------------------------------------------

    private void validateCorrectness(final Node root, final Fixture fixture, final Scratch s) {
        final double rotatedLogL = computeFullGradient(root, fixture, s);
        final double referenceLogL = computeReferenceLogLikelihoodNatural(root, fixture, s.snapshot);
        final double logAbsDetT = -s.snapshot.getLogAbsDetR();
        final double expectedRotatedLogL = referenceLogL - TIP_COUNT * logAbsDetT;
        assertRelativelyClose("p=" + fixture.dim + ": rotated logL vs Jacobian-corrected reference logL",
                expectedRotatedLogL, rotatedLogL, RELATIVE_TOL);

        if (fixture.num2x2Blocks > 0) {
            final double analyticRho = s.nativeGradAcc[0];
            final double originalRho = fixture.rhoParam.getParameterValue(0);
            fixture.rhoParam.setParameterValue(0, originalRho + FD_DELTA);
            final double plusRho = computeReferenceLogLikelihoodNatural(root, fixture, s.snapshot);
            fixture.rhoParam.setParameterValue(0, originalRho - FD_DELTA);
            final double minusRho = computeReferenceLogLikelihoodNatural(root, fixture, s.snapshot);
            fixture.rhoParam.setParameterValue(0, originalRho);
            final double numericRho = (plusRho - minusRho) / (2.0 * FD_DELTA);
            assertRelativelyClose("p=" + fixture.dim + ": dLogL/dRho[block0]", numericRho, analyticRho, RELATIVE_TOL);
        }

        final double analyticR01 = s.rGradAcc.get(0, 1);
        final double originalR01 = fixture.rotationParam.getParameterValue(0, 1);
        fixture.rotationParam.setParameterValue(0, 1, originalR01 + FD_DELTA);
        final double plusR = computeReferenceLogLikelihoodNatural(root, fixture, s.snapshot);
        fixture.rotationParam.setParameterValue(0, 1, originalR01 - FD_DELTA);
        final double minusR = computeReferenceLogLikelihoodNatural(root, fixture, s.snapshot);
        fixture.rotationParam.setParameterValue(0, 1, originalR01);
        final double numericR01 = (plusR - minusR) / (2.0 * FD_DELTA);
        assertRelativelyClose("p=" + fixture.dim + ": dLogL/dR[0,1]", numericR01, analyticR01, RELATIVE_TOL);

        final double analyticQ00 = s.qGradFlatAcc[0];
        final double originalQ00 = fixture.diffusion.getParameterValue(0, 0);
        fixture.diffusion.setParameterValue(0, 0, originalQ00 + FD_DELTA);
        final double plusQ = computeReferenceLogLikelihoodNatural(root, fixture, s.snapshot);
        fixture.diffusion.setParameterValue(0, 0, originalQ00 - FD_DELTA);
        final double minusQ = computeReferenceLogLikelihoodNatural(root, fixture, s.snapshot);
        fixture.diffusion.setParameterValue(0, 0, originalQ00);
        final double numericQ00 = (plusQ - minusQ) / (2.0 * FD_DELTA);
        assertRelativelyClose("p=" + fixture.dim + ": dLogL/dQ[0,0]", numericQ00, analyticQ00, RELATIVE_TOL);
    }

    private static void assertRelativelyClose(final String label, final double expected, final double actual, final double tol) {
        final double scale = Math.max(1.0, Math.abs(expected));
        final double relativeError = Math.abs(actual - expected) / scale;
        if (relativeError > tol) {
            throw new AssertionError(label + ": relative error " + relativeError + " exceeds " + tol
                    + " (expected=" + expected + ", actual=" + actual + ")");
        }
    }

    /** Independent dense/natural-basis reference, generalized to an arbitrary balanced tree. */
    private double computeReferenceLogLikelihoodNatural(final Node root, final Fixture fixture, final BlockBasisGlobalSnapshot snapshot) {
        snapshot.ensureCurrent();
        final int dim = fixture.dim;
        final BlockDiagonalSelectionMatrixParameterization production =
                new BlockDiagonalSelectionMatrixParameterization(fixture.blockParameter);
        final CanonicalGaussianMessageOps.Workspace ws = new CanonicalGaussianMessageOps.Workspace(dim);
        final Map<Node, CanonicalGaussianState> pushedUp = new HashMap<>();
        final Map<Node, CanonicalGaussianState> selfMsg = new HashMap<>();

        naturalPostOrder(root, true, fixture, production, ws, pushedUp, selfMsg);

        final DenseMatrix64F vTilde = snapshot.getStationaryCovarianceDBasis();
        final DenseMatrix64F r = snapshot.getRotation();
        final DenseMatrix64F temp = new DenseMatrix64F(dim, dim);
        final DenseMatrix64F vNatural = new DenseMatrix64F(dim, dim);
        CommonOps.mult(r, vTilde, temp);
        CommonOps.multTransB(temp, r, vNatural);

        final double[] priorCovFlat = new double[dim * dim];
        for (int i = 0; i < dim * dim; ++i) {
            priorCovFlat[i] = vNatural.data[i] / KAPPA0;
        }
        final GaussianFormConverter.Workspace gfcWs = new GaussianFormConverter.Workspace();
        gfcWs.ensureDim(dim);
        final CanonicalGaussianState priorState = new CanonicalGaussianState(dim);
        GaussianFormConverter.fillStateFromMoments(fixture.beta, priorCovFlat, dim, gfcWs, priorState);

        final CanonicalGaussianState combined = new CanonicalGaussianState(dim);
        CanonicalGaussianMessageOps.combineStates(selfMsg.get(root), priorState, combined);
        return CanonicalGaussianMessageOps.normalizationShift(combined, ws);
    }

    private void naturalPostOrder(final Node node,
                                  final boolean isRoot,
                                  final Fixture fixture,
                                  final BlockDiagonalSelectionMatrixParameterization production,
                                  final CanonicalGaussianMessageOps.Workspace ws,
                                  final Map<Node, CanonicalGaussianState> pushedUp,
                                  final Map<Node, CanonicalGaussianState> selfMsg) {
        final int dim = fixture.dim;
        if (node.leaf) {
            final CanonicalGaussianTransition transition = new CanonicalGaussianTransition(dim);
            production.fillCanonicalTransition(fixture.diffusion, fixture.beta, node.branchLength, transition);
            final CanonicalGaussianState out = new CanonicalGaussianState(dim);
            CanonicalGaussianMessageOps.conditionOnObservedSecondBlock(transition, fixture.tipValuesNatural[node.tipIndex], out);
            pushedUp.put(node, out);
            return;
        }
        naturalPostOrder(node.left, false, fixture, production, ws, pushedUp, selfMsg);
        naturalPostOrder(node.right, false, fixture, production, ws, pushedUp, selfMsg);
        final CanonicalGaussianState self = new CanonicalGaussianState(dim);
        CanonicalGaussianMessageOps.combineStates(pushedUp.get(node.left), pushedUp.get(node.right), self);
        selfMsg.put(node, self);
        if (!isRoot) {
            final CanonicalGaussianTransition transition = new CanonicalGaussianTransition(dim);
            production.fillCanonicalTransition(fixture.diffusion, fixture.beta, node.branchLength, transition);
            final CanonicalGaussianState out = new CanonicalGaussianState(dim);
            CanonicalGaussianMessageOps.pushBackward(self, transition, ws, out);
            pushedUp.put(node, out);
        }
    }

    private void computePostOrder(final Node node, final Scratch s) {
        if (node.leaf) {
            s.builder.fillCanonicalTransition(s.snapshot, node.branchLength, node.transition);
            CanonicalGaussianMessageOps.conditionOnObservedSecondBlock(
                    node.transition, s.tipValuesTilde[node.tipIndex], node.pushedUpMessage);
            return;
        }
        computePostOrder(node.left, s);
        computePostOrder(node.right, s);
        CanonicalGaussianMessageOps.combineStates(node.left.pushedUpMessage, node.right.pushedUpMessage, node.selfMessage);
        if (node.transition != null) {
            s.builder.fillCanonicalTransition(s.snapshot, node.branchLength, node.transition);
            CanonicalGaussianMessageOps.pushBackward(node.selfMessage, node.transition, s.ws, node.pushedUpMessage);
        }
    }

    private void computeBackwardAndPullback(final Node node, final CanonicalGaussianState outsideForNode, final Scratch s) {
        final Node left = node.left;
        final Node right = node.right;

        CanonicalGaussianMessageOps.combineStates(outsideForNode, right.pushedUpMessage, s.excludingRight);
        CanonicalGaussianMessageOps.combineStates(outsideForNode, left.pushedUpMessage, s.excludingLeft);

        if (left.leaf) {
            accumulateTipPullback(node.nodeFullPosterior, s.tipValuesTilde[left.tipIndex], left.transition, left.branchLength, s);
        } else {
            CanonicalGaussianMessageOps.buildPairPosterior(s.excludingRight, left.transition, left.selfMessage, s.pairState);
            CanonicalBranchMessageContributionUtils.fillFromPairState(s.pairState, s.pairWs, s.contribution);
            CanonicalTransitionAdjointUtils.fillFromCanonicalTransition(left.transition, s.contribution, s.transitionAdjointWs, s.adjoints);
            s.adjointBuilder.accumulateBranchPullback(s.snapshot, left.branchLength, s.adjoints, s.compressedDAcc, s.gSAcc);
            CanonicalGaussianMessageOps.pushForward(s.excludingRight, left.transition, s.ws, left.outsideMessage);
            CanonicalGaussianMessageOps.combineStates(left.selfMessage, left.outsideMessage, left.nodeFullPosterior);
        }

        if (right.leaf) {
            accumulateTipPullback(node.nodeFullPosterior, s.tipValuesTilde[right.tipIndex], right.transition, right.branchLength, s);
        } else {
            CanonicalGaussianMessageOps.buildPairPosterior(s.excludingLeft, right.transition, right.selfMessage, s.pairState);
            CanonicalBranchMessageContributionUtils.fillFromPairState(s.pairState, s.pairWs, s.contribution);
            CanonicalTransitionAdjointUtils.fillFromCanonicalTransition(right.transition, s.contribution, s.transitionAdjointWs, s.adjoints);
            s.adjointBuilder.accumulateBranchPullback(s.snapshot, right.branchLength, s.adjoints, s.compressedDAcc, s.gSAcc);
            CanonicalGaussianMessageOps.pushForward(s.excludingLeft, right.transition, s.ws, right.outsideMessage);
            CanonicalGaussianMessageOps.combineStates(right.selfMessage, right.outsideMessage, right.nodeFullPosterior);
        }

        if (!left.leaf) {
            computeBackwardAndPullback(left, left.outsideMessage, s);
        }
        if (!right.leaf) {
            computeBackwardAndPullback(right, right.outsideMessage, s);
        }
    }

    /**
     * Fisher contribution for a tip with degenerate (observed, zero-variance) child
     * value, plus the tip-transform's own R-gradient contribution — mirrors {@code
     * BlockBasisCanonicalAdjointBuilderTest#accumulateTipPullback}.
     */
    private void accumulateTipPullback(final CanonicalGaussianState parentPosterior,
                                       final double[] yObserved,
                                       final CanonicalGaussianTransition transition,
                                       final double dt,
                                       final Scratch s) {
        final int dim = s.dim;
        GaussianFormConverter.fillMomentsFromState(parentPosterior, s.parentMean, s.parentCovFlat, dim, s.gfcWs);

        for (int i = 0; i < dim; ++i) {
            s.contribution.dLogL_dInformationX[i] = s.parentMean[i];
            s.contribution.dLogL_dInformationY[i] = yObserved[i];
            for (int j = 0; j < dim; ++j) {
                final double xxMoment = s.parentCovFlat[i * dim + j] + s.parentMean[i] * s.parentMean[j];
                s.contribution.dLogL_dPrecisionXX[i * dim + j] = -0.5 * xxMoment;
                s.contribution.dLogL_dPrecisionXY[i * dim + j] = -0.5 * s.parentMean[i] * yObserved[j];
                s.contribution.dLogL_dPrecisionYX[i * dim + j] = -0.5 * yObserved[i] * s.parentMean[j];
                s.contribution.dLogL_dPrecisionYY[i * dim + j] = -0.5 * yObserved[i] * yObserved[j];
            }
        }
        s.contribution.dLogL_dLogNormalizer = -1.0;

        for (int k = 0; k < dim; ++k) {
            double precisionYYy = 0.0;
            double precisionXYTransposeMean = 0.0;
            for (int i = 0; i < dim; ++i) {
                precisionYYy += transition.precisionYY[k * dim + i] * yObserved[i];
                precisionXYTransposeMean += transition.precisionXY[i * dim + k] * s.parentMean[i];
            }
            s.dLogLdY[k] = -precisionXYTransposeMean - precisionYYy + transition.informationY[k];
        }
        s.adjointBuilder.addTipTransformRGradient(s.snapshot, yObserved, s.dLogLdY, s.rGradAcc);

        CanonicalTransitionAdjointUtils.fillFromCanonicalTransition(transition, s.contribution, s.transitionAdjointWs, s.adjoints);
        s.adjointBuilder.accumulateBranchPullback(s.snapshot, dt, s.adjoints, s.compressedDAcc, s.gSAcc);
    }

    /** The root prior's own direct dependence on {@code Vtilde} (mean is zero regardless of Vtilde). */
    private void accumulateRootPriorVtildePullback(final CanonicalGaussianState rootPosterior,
                                                    final CanonicalGaussianState priorState,
                                                    final Scratch s) {
        final int dim = s.dim;
        GaussianFormConverter.fillMomentsFromState(rootPosterior, s.parentMean, s.parentCovFlat, dim, s.gfcWs);

        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                s.dLogLdPriorCov[i * dim + j] = -0.5 * (s.parentCovFlat[i * dim + j] + s.parentMean[i] * s.parentMean[j]);
            }
        }

        System.arraycopy(priorState.precision, 0, s.priorPrecisionWrap.data, 0, dim * dim);
        System.arraycopy(s.dLogLdPriorCov, 0, s.dLogLdPriorCovWrap.data, 0, dim * dim);
        CommonOps.mult(s.priorPrecisionWrap, s.dLogLdPriorCovWrap, s.temp1);
        CommonOps.mult(s.temp1, s.priorPrecisionWrap, s.temp2);

        for (int i = 0; i < dim * dim; ++i) {
            s.rootPriorResult[i] = -(1.0 / KAPPA0) * (s.temp2.data[i] + 0.5 * s.priorPrecisionWrap.data[i]);
        }
        s.adjointBuilder.accumulateRootPriorPullback(s.rootPriorResult, s.gSAcc);
    }

    // ------------------------------------------------------------------
    // Per-iteration parameter perturbation (forces ensureCurrent() to redo its work).
    // ------------------------------------------------------------------

    private void perturbParameters(final int iteration, final Fixture fixture) {
        final double eps = 1.0e-7 * (iteration + 1);
        for (int b = 0; b < fixture.num2x2Blocks; ++b) {
            fixture.rhoParam.setParameterValue(b, fixture.baseRho[b] + eps);
            fixture.thetaParam.setParameterValue(b, fixture.baseTheta[b] + eps);
            fixture.tParam.setParameterValue(b, fixture.baseT[b] + eps);
        }
        if (fixture.scalarParam.getDimension() > 0) {
            fixture.scalarParam.setParameterValue(0, fixture.baseScalar[0] + eps);
        }
        final int dim = fixture.dim;
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                if (i != j) {
                    fixture.rotationParam.setParameterValue(i, j, fixture.baseRotation[i * dim + j] + eps);
                }
            }
        }
    }

    // ------------------------------------------------------------------
    // Fixture / tree construction.
    // ------------------------------------------------------------------

    private static final class Fixture {
        final int dim;
        final int num2x2Blocks;
        final MatrixParameter rotationParam;
        final double[] baseRotation;
        final Parameter scalarParam;
        final double[] baseScalar;
        final Parameter rhoParam;
        final double[] baseRho;
        final Parameter thetaParam;
        final double[] baseTheta;
        final Parameter tParam;
        final double[] baseT;
        final BlockDiagonalPolarStableMatrixParameter blockParameter;
        final MatrixParameter diffusion;
        final double[] beta;
        final double[][] tipValuesNatural;

        Fixture(final int dim, final int num2x2Blocks,
               final MatrixParameter rotationParam, final double[] baseRotation,
               final Parameter scalarParam, final double[] baseScalar,
               final Parameter rhoParam, final double[] baseRho,
               final Parameter thetaParam, final double[] baseTheta,
               final Parameter tParam, final double[] baseT,
               final BlockDiagonalPolarStableMatrixParameter blockParameter,
               final MatrixParameter diffusion, final double[] beta, final double[][] tipValuesNatural) {
            this.dim = dim;
            this.num2x2Blocks = num2x2Blocks;
            this.rotationParam = rotationParam;
            this.baseRotation = baseRotation;
            this.scalarParam = scalarParam;
            this.baseScalar = baseScalar;
            this.rhoParam = rhoParam;
            this.baseRho = baseRho;
            this.thetaParam = thetaParam;
            this.baseTheta = baseTheta;
            this.tParam = tParam;
            this.baseT = baseT;
            this.blockParameter = blockParameter;
            this.diffusion = diffusion;
            this.beta = beta;
            this.tipValuesNatural = tipValuesNatural;
        }
    }

    private Fixture buildFixture(final int dim, final Random random) {
        final boolean oddDim = (dim & 1) == 1;
        final int num2x2Blocks = dim / 2;

        final double[] baseRotation = new double[dim * dim];
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                baseRotation[i * dim + j] = (i == j) ? 1.0 : (random.nextDouble() - 0.5) * (0.3 / dim);
            }
        }
        final MatrixParameter rotationParam = new MatrixParameter("R", dim, dim);
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                rotationParam.setParameterValue(i, j, baseRotation[i * dim + j]);
            }
        }

        final double[] baseRho = new double[num2x2Blocks];
        final double[] baseTheta = new double[num2x2Blocks];
        final double[] baseT = new double[num2x2Blocks];
        for (int b = 0; b < num2x2Blocks; ++b) {
            baseRho[b] = 0.5 + random.nextDouble();
            baseTheta[b] = (random.nextDouble() - 0.5) * 1.2;
            baseT[b] = (random.nextDouble() - 0.5) * 0.4;
        }
        final Parameter rhoParam = new Parameter.Default(baseRho.clone());
        final Parameter thetaParam = new Parameter.Default(baseTheta.clone());
        final Parameter tParam = new Parameter.Default(baseT.clone());

        final double[] baseScalar = oddDim ? new double[]{0.5 + random.nextDouble()} : new double[0];
        final Parameter scalarParam = new Parameter.Default(baseScalar.length);
        if (oddDim) {
            scalarParam.setParameterValue(0, baseScalar[0]);
        }

        final BlockDiagonalPolarStableMatrixParameter blockParameter = new BlockDiagonalPolarStableMatrixParameter(
                "benchmarkBlock" + dim, rotationParam, scalarParam, rhoParam, thetaParam, tParam);

        final MatrixParameter diffusion = randomSpdMatrix(dim, random);

        final double[] beta = new double[dim];
        for (int i = 0; i < dim; ++i) {
            beta[i] = random.nextGaussian();
        }

        final double[][] tipValuesNatural = new double[TIP_COUNT][dim];
        for (int t = 0; t < TIP_COUNT; ++t) {
            for (int i = 0; i < dim; ++i) {
                tipValuesNatural[t][i] = 3.0 * Math.sin(0.13 * (t + 1) * (i + 1));
            }
        }

        return new Fixture(dim, num2x2Blocks, rotationParam, baseRotation, scalarParam, baseScalar,
                rhoParam, baseRho, thetaParam, baseTheta, tParam, baseT, blockParameter, diffusion, beta, tipValuesNatural);
    }

    private static MatrixParameter randomSpdMatrix(final int dim, final Random random) {
        final double[][] a = new double[dim][dim];
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                a[i][j] = random.nextDouble() - 0.5;
            }
        }
        final MatrixParameter q = new MatrixParameter("Q", dim, dim);
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                double sum = 0.0;
                for (int k = 0; k < dim; ++k) {
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

    private static void transformToTilde(final BlockBasisGlobalSnapshot snapshot,
                                         final double[] y,
                                         final double[] beta,
                                         final int dim,
                                         final double[] out) {
        final DenseMatrix64F t = snapshot.getRotationInverse();
        for (int i = 0; i < dim; ++i) {
            double sum = 0.0;
            for (int j = 0; j < dim; ++j) {
                sum += t.get(i, j) * (y[j] - beta[j]);
            }
            out[i] = sum;
        }
    }

    private static final class Node {
        Node left;
        Node right;
        boolean leaf;
        int tipIndex;
        double branchLength;

        CanonicalGaussianTransition transition;
        CanonicalGaussianState pushedUpMessage;
        CanonicalGaussianState selfMessage;
        CanonicalGaussianState outsideMessage;
        CanonicalGaussianState nodeFullPosterior;
    }

    private Node buildBalancedTree(final Random random) {
        return buildSubtree(0, TIP_COUNT, random, true);
    }

    private Node buildSubtree(final int start, final int end, final Random random, final boolean isRoot) {
        final Node node = new Node();
        if (end - start == 1) {
            node.leaf = true;
            node.tipIndex = start;
        } else {
            node.leaf = false;
            final int mid = start + (end - start) / 2;
            node.left = buildSubtree(start, mid, random, false);
            node.right = buildSubtree(mid, end, random, false);
        }
        if (!isRoot) {
            node.branchLength = 0.05 + random.nextDouble() * 0.55;
        }
        return node;
    }

    /** Allocates every per-node state object once (sized to {@code dim}), reused across iterations. */
    private void allocateNodeStates(final Node node, final int dim, final boolean isRoot) {
        if (!isRoot) {
            node.transition = new CanonicalGaussianTransition(dim);
            node.pushedUpMessage = new CanonicalGaussianState(dim);
        }
        if (!node.leaf) {
            node.selfMessage = new CanonicalGaussianState(dim);
            node.nodeFullPosterior = new CanonicalGaussianState(dim);
            if (!isRoot) {
                node.outsideMessage = new CanonicalGaussianState(dim);
            }
            allocateNodeStates(node.left, dim, false);
            allocateNodeStates(node.right, dim, false);
        }
    }

    private static boolean isFinite(final double[] values) {
        for (final double value : values) {
            if (!Double.isFinite(value)) {
                return false;
            }
        }
        return true;
    }

    // ------------------------------------------------------------------
    // Reused scratch (allocated once per dimension, outside the timed loop).
    // ------------------------------------------------------------------

    private final class Scratch {
        final int dim;
        final BlockBasisGlobalSnapshot snapshot;
        final BlockBasisCanonicalTransitionBuilder builder;
        final BlockBasisCanonicalAdjointBuilder adjointBuilder;
        final BlockDiagonalPolarStableMatrixParameter blockParameter;

        final double[][] tipValuesTilde;
        final double[] zeroMean;
        final double[] priorCovFlat;
        final CanonicalGaussianState priorState;

        final double[] compressedDAcc;
        final DenseMatrix64F gSAcc;
        final DenseMatrix64F rGradAcc;
        final double[] qGradFlatAcc;
        final double[] nativeGradAcc;

        final CanonicalGaussianMessageOps.Workspace ws;
        final GaussianFormConverter.Workspace gfcWs;
        final CanonicalTransitionAdjointUtils.Workspace transitionAdjointWs;
        final CanonicalLocalTransitionAdjoints adjoints;
        final CanonicalBranchMessageContributionUtils.Workspace pairWs;
        final CanonicalBranchMessageContribution contribution;
        final CanonicalGaussianState pairState;
        final CanonicalGaussianState excludingRight;
        final CanonicalGaussianState excludingLeft;

        final double[] parentMean;
        final double[] parentCovFlat;
        final double[] dLogLdY;
        final double[] dLogLdPriorCov;
        final double[] rootPriorResult;
        final DenseMatrix64F priorPrecisionWrap;
        final DenseMatrix64F dLogLdPriorCovWrap;
        final DenseMatrix64F temp1;
        final DenseMatrix64F temp2;

        Scratch(final int dim, final Fixture fixture) {
            this.dim = dim;
            this.blockParameter = fixture.blockParameter;
            this.snapshot = new BlockBasisGlobalSnapshot(fixture.blockParameter, fixture.diffusion);
            this.snapshot.ensureCurrent();
            this.builder = new BlockBasisCanonicalTransitionBuilder(dim, snapshot.getBlockStarts(), snapshot.getBlockSizes());
            this.adjointBuilder = new BlockBasisCanonicalAdjointBuilder(dim, snapshot.getBlockStarts(), snapshot.getBlockSizes());

            this.tipValuesTilde = new double[TIP_COUNT][dim];
            this.zeroMean = new double[dim];
            this.priorCovFlat = new double[dim * dim];
            this.priorState = new CanonicalGaussianState(dim);

            this.compressedDAcc = new double[fixture.blockParameter.getCompressedDDimension()];
            this.gSAcc = new DenseMatrix64F(dim, dim);
            this.rGradAcc = new DenseMatrix64F(dim, dim);
            this.qGradFlatAcc = new double[dim * dim];
            this.nativeGradAcc = new double[fixture.blockParameter.getBlockDiagonalNParameters()];

            this.ws = new CanonicalGaussianMessageOps.Workspace(dim);
            this.gfcWs = new GaussianFormConverter.Workspace();
            this.gfcWs.ensureDim(dim);
            this.transitionAdjointWs = new CanonicalTransitionAdjointUtils.Workspace(dim);
            this.adjoints = new CanonicalLocalTransitionAdjoints(dim);
            this.pairWs = new CanonicalBranchMessageContributionUtils.Workspace(dim);
            this.contribution = new CanonicalBranchMessageContribution(dim);
            this.pairState = new CanonicalGaussianState(2 * dim);
            this.excludingRight = new CanonicalGaussianState(dim);
            this.excludingLeft = new CanonicalGaussianState(dim);

            this.parentMean = new double[dim];
            this.parentCovFlat = new double[dim * dim];
            this.dLogLdY = new double[dim];
            this.dLogLdPriorCov = new double[dim * dim];
            this.rootPriorResult = new double[dim * dim];
            this.priorPrecisionWrap = new DenseMatrix64F(dim, dim);
            this.dLogLdPriorCovWrap = new DenseMatrix64F(dim, dim);
            this.temp1 = new DenseMatrix64F(dim, dim);
            this.temp2 = new DenseMatrix64F(dim, dim);
        }
    }

}

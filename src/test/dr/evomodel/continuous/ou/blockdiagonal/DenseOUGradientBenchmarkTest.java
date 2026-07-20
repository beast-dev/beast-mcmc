package test.dr.evomodel.continuous.ou.blockdiagonal;

import dr.evolution.io.Importer;
import dr.evolution.io.NewickImporter;
import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.DefaultBranchRateModel;
import dr.evomodel.continuous.MultivariateDiffusionModel;
import dr.evomodel.continuous.MultivariateElasticModel;
import dr.evomodel.tree.DefaultTreeModel;
import dr.evomodel.treedatalikelihood.continuous.ConjugateRootTraitPrior;
import dr.evomodel.treedatalikelihood.continuous.canonical.CanonicalRootPrior;
import dr.evomodel.treedatalikelihood.continuous.canonical.SequentialCanonicalOUMessagePasser;
import dr.evomodel.treedatalikelihood.continuous.canonical.adapter.CanonicalConjugateRootPriorAdapter;
import dr.evomodel.treedatalikelihood.continuous.canonical.adapter.HomogeneousCanonicalOUBranchTransitionProvider;
import dr.evomodel.treedatalikelihood.continuous.observationmodel.CanonicalTipObservation;
import dr.inference.model.BlockDiagonalPolarStableMatrixParameter;
import dr.inference.model.MatrixParameter;
import dr.inference.model.Parameter;
import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

import java.io.IOException;
import java.util.Random;

/**
 * Rotated/block-basis OU work, gradient-benchmark counterpart (see
 * private-code/rotatedOU/rotated_ou_block_basis_plan.md and
 * private-code/rotatedOU/phase3_gradient_benchmark.md).
 *
 * <p>Standalone (not part of {@code ant junit-canonical}) wall-clock benchmark of one
 * <b>full gradient evaluation</b> under the <em>existing, unmodified, dense-rotation</em>
 * canonical OU backend ({@code BlockDiagonalSelectionMatrixParameterization}, reached via
 * {@code MultivariateElasticModel} + {@code HomogeneousCanonicalOUBranchTransitionProvider}
 * + {@code SequentialCanonicalOUMessagePasser}) — i.e. the "previous approach" that does
 * <b>not</b> rotate tip values into a {@code Ytilde}-basis, working entirely in natural
 * coordinates with a per-branch dense {@code S = R*B*R^{-1}} sandwich. Same balanced 50-tip
 * tree and {@code p = 4, 8, 16} as {@link BlockBasisGradientBenchmarkTest}, for direct,
 * like-for-like comparison.
 *
 * <p>Unlike the new backend's benchmark, this one needs no hand-rolled tree recursion: the
 * dense backend already has a generalized, production tree-traversal/gradient path
 * ({@code SequentialCanonicalOUMessagePasser.computePostOrderLogLikelihood/computePreOrder/
 * computeJointGradients}), the same one {@code RotatedOUBaselineBenchmarkTest} already uses
 * for its own 128/1024-tip synthetic cases.
 */
public final class DenseOUGradientBenchmarkTest extends TestCase {

    private static final int TIP_COUNT = 100;
    private static final int[] DIMENSIONS = {4, 8, 16};
    private static final int WARMUP_ITERS = 200;
    private static final int MEASURE_ITERS = 1000;
    private static final double KAPPA0 = 10.0;
    private static final long SEED = 42L;
    private static final double FD_DELTA = 1.0e-6;
    private static final double RELATIVE_TOL = 5.0e-4;

    public DenseOUGradientBenchmarkTest(final String name) {
        super(name);
    }

    public static Test suite() {
        return new TestSuite(DenseOUGradientBenchmarkTest.class);
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
        final Random treeRandom = new Random(SEED + dim);
        final DefaultTreeModel treeModel = buildBalancedTreeModel(dim, treeRandom);

        final Random fixtureRandom = new Random(SEED + dim + 1);
        final Fixture fixture = buildFixture(dim, fixtureRandom);

        final HomogeneousCanonicalOUBranchTransitionProvider provider = new HomogeneousCanonicalOUBranchTransitionProvider(
                treeModel, new MultivariateElasticModel(fixture.blockSelection), fixture.diffusionModel,
                fixture.muParam, new DefaultBranchRateModel());

        final SequentialCanonicalOUMessagePasser passer = new SequentialCanonicalOUMessagePasser(treeModel, dim);
        for (int t = 0; t < TIP_COUNT; ++t) {
            final CanonicalTipObservation observation = new CanonicalTipObservation(dim);
            observation.setObserved(fixture.tipValues[t]);
            passer.setTipObservation(t, observation);
        }

        final ConjugateRootTraitPrior rootPriorRaw = new ConjugateRootTraitPrior(
                new Parameter.Default(fixture.beta.clone()), new Parameter.Default(KAPPA0));
        final CanonicalRootPrior rootPrior = new CanonicalConjugateRootPriorAdapter(rootPriorRaw, dim);

        final int nativeBlockDim = fixture.blockSelection.getBlockDiagonalNParameters();
        final int gradADim = nativeBlockDim + dim * dim;
        final double[] gradA = new double[gradADim];
        final double[] gradQ = new double[dim * dim];
        final double[] gradMu = new double[dim];

        validateCorrectness(passer, provider, rootPrior, fixture, gradA, gradQ, gradMu);

        double logL = Double.NaN;
        for (int i = 0; i < WARMUP_ITERS; ++i) {
            perturbParameters(i, fixture);
            logL = passer.computePostOrderLogLikelihood(provider, rootPrior);
            passer.computePreOrder(provider, rootPrior);
            passer.computeJointGradients(provider, gradA, gradQ, gradMu);
        }
        if (!Double.isFinite(logL)) {
            throw new AssertionError("Non-finite benchmark log-likelihood at p=" + dim + ": " + logL);
        }

        final long start = System.nanoTime();
        for (int i = 0; i < MEASURE_ITERS; ++i) {
            perturbParameters(WARMUP_ITERS + i, fixture);
            logL = passer.computePostOrderLogLikelihood(provider, rootPrior);
            passer.computePreOrder(provider, rootPrior);
            passer.computeJointGradients(provider, gradA, gradQ, gradMu);
        }
        final long elapsedNs = System.nanoTime() - start;
        final double nsPerGradient = elapsedNs / (double) MEASURE_ITERS;

        if (!isFinite(gradA) || !isFinite(gradQ) || !isFinite(gradMu)) {
            throw new AssertionError("Non-finite benchmark gradient at p=" + dim);
        }

        final int edgeCount = 2 * TIP_COUNT - 2;
        final StringBuilder sb = new StringBuilder();
        sb.append("Case p=").append(dim).append(", ").append(TIP_COUNT).append(" tips, ")
                .append(edgeCount).append(" edges (dense/previous approach, tips not rotated)\n");
        sb.append(String.format("  logL: %.6f%n", logL));
        sb.append(String.format("  full gradient (post-order + pre-order + joint gradient): %.1f ns (%.3f us)%n",
                nsPerGradient, nsPerGradient / 1000.0));
        sb.append(String.format("  per edge: %.1f ns%n", nsPerGradient / edgeCount));
        sb.append(String.format("  throughput: %.1f gradients/sec%n", 1.0e9 / nsPerGradient));
        sb.append('\n');
        return sb.toString();
    }

    // ------------------------------------------------------------------
    // Correctness validation (FD spot checks against this same setup's own likelihood --
    // no cross-backend reference needed, since computePostOrderLogLikelihood/
    // computeJointGradients are the existing, already-tested production path).
    // ------------------------------------------------------------------

    private void validateCorrectness(final SequentialCanonicalOUMessagePasser passer,
                                     final HomogeneousCanonicalOUBranchTransitionProvider provider,
                                     final CanonicalRootPrior rootPrior,
                                     final Fixture fixture,
                                     final double[] gradA,
                                     final double[] gradQ,
                                     final double[] gradMu) {
        passer.computePostOrderLogLikelihood(provider, rootPrior);
        passer.computePreOrder(provider, rootPrior);
        passer.computeJointGradients(provider, gradA, gradQ, gradMu);

        if (fixture.num2x2Blocks > 0) {
            final double analyticRho = gradA[0];
            final double originalRho = fixture.rhoParam.getParameterValue(0);
            fixture.rhoParam.setParameterValue(0, originalRho + FD_DELTA);
            final double plusRho = passer.computePostOrderLogLikelihood(provider, rootPrior);
            fixture.rhoParam.setParameterValue(0, originalRho - FD_DELTA);
            final double minusRho = passer.computePostOrderLogLikelihood(provider, rootPrior);
            fixture.rhoParam.setParameterValue(0, originalRho);
            final double numericRho = (plusRho - minusRho) / (2.0 * FD_DELTA);
            assertRelativelyClose("p=" + fixture.dim + ": dLogL/dRho[block0] (dense)", numericRho, analyticRho, RELATIVE_TOL);
        }

        final int nativeBlockDim = fixture.blockSelection.getBlockDiagonalNParameters();
        // gradA's rotation tail is flattened column-major: entry (row,col) lives at col*dim+row.
        final double analyticR01 = gradA[nativeBlockDim + 1 * fixture.dim + 0];
        final double originalR01 = fixture.rotationParam.getParameterValue(0, 1);
        fixture.rotationParam.setParameterValue(0, 1, originalR01 + FD_DELTA);
        final double plusR = passer.computePostOrderLogLikelihood(provider, rootPrior);
        fixture.rotationParam.setParameterValue(0, 1, originalR01 - FD_DELTA);
        final double minusR = passer.computePostOrderLogLikelihood(provider, rootPrior);
        fixture.rotationParam.setParameterValue(0, 1, originalR01);
        final double numericR01 = (plusR - minusR) / (2.0 * FD_DELTA);
        assertRelativelyClose("p=" + fixture.dim + ": dLogL/dR[0,1] (dense, col-major)", numericR01, analyticR01, RELATIVE_TOL);

        final double analyticQ00 = gradQ[0];
        final double originalQ00 = fixture.qMatrix.getParameterValue(0, 0);
        fixture.qMatrix.setParameterValue(0, 0, originalQ00 + FD_DELTA);
        syncPrecisionFromVariance(fixture.qMatrix, fixture.precisionMatrix);
        final double plusQ = passer.computePostOrderLogLikelihood(provider, rootPrior);
        fixture.qMatrix.setParameterValue(0, 0, originalQ00 - FD_DELTA);
        syncPrecisionFromVariance(fixture.qMatrix, fixture.precisionMatrix);
        final double minusQ = passer.computePostOrderLogLikelihood(provider, rootPrior);
        fixture.qMatrix.setParameterValue(0, 0, originalQ00);
        syncPrecisionFromVariance(fixture.qMatrix, fixture.precisionMatrix);
        final double numericQ00 = (plusQ - minusQ) / (2.0 * FD_DELTA);
        assertRelativelyClose("p=" + fixture.dim + ": dLogL/dQ[0,0] (dense)", numericQ00, analyticQ00, RELATIVE_TOL);

        passer.computePostOrderLogLikelihood(provider, rootPrior);
        passer.computePreOrder(provider, rootPrior);
        passer.computeJointGradients(provider, gradA, gradQ, gradMu);
    }

    private static void assertRelativelyClose(final String label, final double expected, final double actual, final double tol) {
        final double scale = Math.max(1.0, Math.abs(expected));
        final double relativeError = Math.abs(actual - expected) / scale;
        if (relativeError > tol) {
            throw new AssertionError(label + ": relative error " + relativeError + " exceeds " + tol
                    + " (expected=" + expected + ", actual=" + actual + ")");
        }
    }

    private static void syncPrecisionFromVariance(final MatrixParameter qMatrix, final MatrixParameter precisionMatrix) {
        final int dim = qMatrix.getRowDimension();
        final DenseMatrix64F variance = new DenseMatrix64F(dim, dim);
        final DenseMatrix64F precision = new DenseMatrix64F(dim, dim);
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                variance.set(i, j, qMatrix.getParameterValue(i, j));
            }
        }
        if (!CommonOps.invert(variance, precision)) {
            throw new IllegalArgumentException("Variance matrix is singular during gradient benchmark.");
        }
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                precisionMatrix.setParameterValueQuietly(i, j, precision.get(i, j));
            }
        }
        precisionMatrix.fireParameterChangedEvent();
    }

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

    private static boolean isFinite(final double[] values) {
        for (final double value : values) {
            if (!Double.isFinite(value)) {
                return false;
            }
        }
        return true;
    }

    // ------------------------------------------------------------------
    // Fixture / tree construction (mirrors BlockBasisGradientBenchmarkTest's Fixture,
    // self-contained here since the dense backend doesn't need the BlockBasis* classes).
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
        final BlockDiagonalPolarStableMatrixParameter blockSelection;
        final MatrixParameter qMatrix;
        final MatrixParameter precisionMatrix;
        final MultivariateDiffusionModel diffusionModel;
        final Parameter muParam;
        final double[] beta;
        final double[][] tipValues;

        Fixture(final int dim, final int num2x2Blocks,
               final MatrixParameter rotationParam, final double[] baseRotation,
               final Parameter scalarParam, final double[] baseScalar,
               final Parameter rhoParam, final double[] baseRho,
               final Parameter thetaParam, final double[] baseTheta,
               final Parameter tParam, final double[] baseT,
               final BlockDiagonalPolarStableMatrixParameter blockSelection,
               final MatrixParameter qMatrix, final MatrixParameter precisionMatrix,
               final MultivariateDiffusionModel diffusionModel,
               final Parameter muParam, final double[] beta, final double[][] tipValues) {
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
            this.blockSelection = blockSelection;
            this.qMatrix = qMatrix;
            this.precisionMatrix = precisionMatrix;
            this.diffusionModel = diffusionModel;
            this.muParam = muParam;
            this.beta = beta;
            this.tipValues = tipValues;
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

        final BlockDiagonalPolarStableMatrixParameter blockSelection = new BlockDiagonalPolarStableMatrixParameter(
                "denseBenchmarkBlock" + dim, rotationParam, scalarParam, rhoParam, thetaParam, tParam);

        final MatrixParameter qMatrix = randomSpdMatrix(dim, random);
        final MatrixParameter precisionMatrix = new MatrixParameter("precision", dim, dim);
        syncPrecisionFromVariance(qMatrix, precisionMatrix);
        final MultivariateDiffusionModel diffusionModel = new MultivariateDiffusionModel(precisionMatrix);

        final double[] beta = new double[dim];
        for (int i = 0; i < dim; ++i) {
            beta[i] = random.nextGaussian();
        }
        final Parameter muParam = new Parameter.Default(beta.clone());

        final double[][] tipValues = new double[TIP_COUNT][dim];
        for (int t = 0; t < TIP_COUNT; ++t) {
            for (int i = 0; i < dim; ++i) {
                tipValues[t][i] = 3.0 * Math.sin(0.13 * (t + 1) * (i + 1));
            }
        }

        return new Fixture(dim, num2x2Blocks, rotationParam, baseRotation, scalarParam, baseScalar,
                rhoParam, baseRho, thetaParam, baseTheta, tParam, baseT, blockSelection,
                qMatrix, precisionMatrix, diffusionModel, muParam, beta, tipValues);
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

    private DefaultTreeModel buildBalancedTreeModel(final int dim, final Random random) {
        try {
            final Tree tree = new NewickImporter(buildBalancedNewick(0, TIP_COUNT, random, true)).importTree(null);
            return new DefaultTreeModel("denseBenchmarkTree" + dim, tree);
        } catch (IOException | Importer.ImportException e) {
            throw new RuntimeException("Failed to build synthetic benchmark tree", e);
        }
    }

    private static String buildBalancedNewick(final int start, final int end, final Random random, final boolean isRoot) {
        final String body;
        if (end - start == 1) {
            body = "taxon" + start;
        } else {
            final int mid = start + (end - start) / 2;
            body = "(" + buildBalancedNewick(start, mid, random, false)
                    + "," + buildBalancedNewick(mid, end, random, false) + ")";
        }
        if (isRoot) {
            return body + ";";
        }
        final double branchLength = 0.05 + random.nextDouble() * 0.55;
        return body + ":" + branchLength;
    }
}

package test.dr.evomodel.continuous.ou.blockdiagonal;

import dr.evolution.io.Importer;
import dr.evolution.io.NewickImporter;
import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.DefaultBranchRateModel;
import dr.evomodel.continuous.MultivariateDiffusionModel;
import dr.evomodel.continuous.MultivariateElasticModel;
import dr.evomodel.continuous.ou.blockdiagonal.BlockBasisCanonicalAdjointBuilder;
import dr.evomodel.continuous.ou.blockdiagonal.BlockBasisCanonicalTransitionBuilder;
import dr.evomodel.continuous.ou.blockdiagonal.BlockBasisGlobalSnapshot;
import dr.evomodel.tree.DefaultTreeModel;
import dr.evomodel.treedatalikelihood.continuous.ConjugateRootTraitPrior;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.CanonicalRootPrior;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.SequentialCanonicalOUMessagePasser;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.adapter.CanonicalConjugateRootPriorAdapter;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.adapter.HomogeneousCanonicalOUBranchTransitionProvider;
import dr.evomodel.continuous.canonical.math.GaussianFormConverter;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.message.CanonicalBranchMessageContribution;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.message.CanonicalBranchMessageContributionUtils;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.message.CanonicalGaussianMessageOps;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.message.CanonicalGaussianState;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.message.CanonicalGaussianTransition;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.message.CanonicalLocalTransitionAdjoints;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.message.CanonicalTransitionAdjointUtils;
import dr.evomodel.treedatalikelihood.continuousNew.observationmodel.CanonicalTipObservation;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.BlockDiagonalPolarStableMatrixParameter;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Likelihood;
import dr.inference.model.MatrixParameter;
import dr.inference.model.Parameter;
import dr.inference.operators.AdaptationMode;
import dr.inference.operators.hmc.HamiltonianMonteCarloOperator;
import dr.inference.operators.hmc.MassPreconditionScheduler;
import dr.inference.operators.hmc.MassPreconditioner;
import dr.inference.operators.hmc.MassPreconditioningOptions;
import dr.math.MathUtils;
import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

import java.io.IOException;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.Locale;
import java.util.Map;
import java.util.Random;

/**
 * Standalone profile for the rotated/block-basis OU prototype and the older
 * dense homogeneous canonical OU path when both are driven by BEAST's vanilla
 * {@link HamiltonianMonteCarloOperator}.
 *
 * <p>The block-basis machinery is not XML/HMC-wired yet, so this test wraps it
 * directly as a {@link GradientWrtParameterProvider}. The dense side uses the
 * production homogeneous canonical OU provider and message passer. Both sides
 * use a matched 100-tip, 16D synthetic fixture by default.</p>
 */
public final class RotatedOUHmcProfileTest extends TestCase {

    private static final int TIP_COUNT = Integer.getInteger("beast.profile.rotatedOU.tips", 100);
    private static final int DIMENSION = Integer.getInteger("beast.profile.rotatedOU.dim", 16);
    private static final int HMC_STEPS = Integer.getInteger("beast.profile.rotatedOU.hmcSteps", 10);
    private static final int WARMUP_OPERATIONS = Integer.getInteger("beast.profile.rotatedOU.warmup", 10);
    private static final int MEASURE_OPERATIONS = Integer.getInteger("beast.profile.rotatedOU.measure", 50);
    private static final double STEP_SIZE = Double.parseDouble(System.getProperty(
            "beast.profile.rotatedOU.stepSize", "1.0e-8"));
    private static final double KAPPA0 = 10.0;
    private static final long FIXTURE_SEED = 42L;
    private static final long HMC_SEED = 666L;

    public RotatedOUHmcProfileTest(final String name) {
        super(name);
    }

    public static Test suite() {
        return new TestSuite(RotatedOUHmcProfileTest.class);
    }

    public static void main(final String[] args) {
        junit.textui.TestRunner.run(suite());
    }

    public void testProfileReportIsFinite() {
        final FixtureData data = buildFixtureData(DIMENSION);

        final DenseProfileEvaluator dense = new DenseProfileEvaluator(data);
        final BlockBasisProfileEvaluator blockBasis = new BlockBasisProfileEvaluator(data);

        runHmcProfile(dense);
        runHmcProfile(blockBasis);

        final String report = buildReport(dense, blockBasis);
        System.out.print(report);

        assertFalse(report.contains("NaN"));
        assertFalse(report.contains("Infinity"));
        assertTrue(dense.stats.gradient.calls > 0);
        assertTrue(blockBasis.stats.gradient.calls > 0);
    }

    private static void runHmcProfile(final BaseEvaluator evaluator) {
        evaluator.resetToBase();
        MathUtils.setSeed(HMC_SEED);
        HamiltonianMonteCarloOperator operator = createOperator(evaluator);
        for (int i = 0; i < WARMUP_OPERATIONS; ++i) {
            runOneHmcOperation(evaluator, operator);
        }

        evaluator.resetToBase();
        evaluator.stats.reset();
        MathUtils.setSeed(HMC_SEED);
        operator = createOperator(evaluator);
        for (int i = 0; i < MEASURE_OPERATIONS; ++i) {
            runOneHmcOperation(evaluator, operator);
        }
    }

    private static void runOneHmcOperation(final BaseEvaluator evaluator,
                                           final HamiltonianMonteCarloOperator operator) {
        final long start = System.nanoTime();
        final double oldLogL = evaluator.getLikelihood().getLogLikelihood();
        final double hastings = operator.operate(evaluator.getLikelihood());
        final double newLogL = evaluator.getLikelihood().getLogLikelihood();
        evaluator.stats.hmcProposal.add(System.nanoTime() - start);
        if (!Double.isFinite(oldLogL) || !Double.isFinite(newLogL) || !Double.isFinite(hastings)) {
            throw new AssertionError(evaluator.label + " produced non-finite HMC state: oldLogL="
                    + oldLogL + " newLogL=" + newLogL + " hastings=" + hastings);
        }
        operator.accept(0.0);
    }

    private static HamiltonianMonteCarloOperator createOperator(final BaseEvaluator evaluator) {
        final MassPreconditioningOptions.Default preconditioningOptions =
                new MassPreconditioningOptions.Default(0, 0, 0, 0, false, null, null);
        final HamiltonianMonteCarloOperator.Options options = new HamiltonianMonteCarloOperator.Options(
                STEP_SIZE,
                HMC_STEPS,
                0.0,
                preconditioningOptions,
                0,
                1.0e-3,
                10,
                0.1,
                0.8,
                HamiltonianMonteCarloOperator.InstabilityHandler.REJECT);
        final MassPreconditioner preconditioner =
                MassPreconditioner.Type.NONE.factory(evaluator, null, options);
        return new HamiltonianMonteCarloOperator(
                AdaptationMode.ADAPTATION_OFF,
                1.0,
                evaluator,
                evaluator.getParameter(),
                null,
                evaluator.mask,
                options,
                preconditioner,
                MassPreconditionScheduler.Type.NONE);
    }

    private static String buildReport(final DenseProfileEvaluator dense,
                                      final BlockBasisProfileEvaluator blockBasis) {
        final StringBuilder sb = new StringBuilder();
        sb.append("Rotated OU HMC profile\n");
        sb.append("  fixture: ").append(TIP_COUNT).append(" tips, p=").append(DIMENSION)
                .append(", edges=").append(2 * TIP_COUNT - 2).append('\n');
        sb.append("  HMC: vanilla, nSteps=").append(HMC_STEPS)
                .append(", stepSize=").append(STEP_SIZE)
                .append(", warmup=").append(WARMUP_OPERATIONS)
                .append(", measured proposals=").append(MEASURE_OPERATIONS)
                .append(", gradientCheckCount=0, preconditioning=none\n\n");
        appendEvaluatorReport(sb, dense);
        appendEvaluatorReport(sb, blockBasis);
        sb.append(String.format(Locale.US,
                "Speed ratios, dense / block-basis: HMC proposal %.3fx, gradient %.3fx, logL %.3fx%n",
                dense.stats.hmcProposal.meanNanos() / blockBasis.stats.hmcProposal.meanNanos(),
                dense.stats.gradient.meanNanos() / blockBasis.stats.gradient.meanNanos(),
                dense.stats.logLikelihood.meanNanos() / blockBasis.stats.logLikelihood.meanNanos()));
        sb.append('\n');
        return sb.toString();
    }

    private static void appendEvaluatorReport(final StringBuilder sb, final BaseEvaluator evaluator) {
        sb.append(evaluator.label).append('\n');
        appendTimer(sb, "HMC proposal total", evaluator.stats.hmcProposal, "proposal");
        appendTimer(sb, "likelihood eval", evaluator.stats.logLikelihood, "call");
        appendTimer(sb, "gradient eval", evaluator.stats.gradient, "call");
        for (final Map.Entry<String, Timer> entry : evaluator.stats.phases.entrySet()) {
            appendTimer(sb, "  " + entry.getKey(), entry.getValue(), "gradient");
        }
        sb.append('\n');
    }

    private static void appendTimer(final StringBuilder sb,
                                    final String label,
                                    final Timer timer,
                                    final String unit) {
        sb.append(String.format(Locale.US,
                "  %-28s %10.3f ms / %s  (calls=%d)%n",
                label + ":",
                timer.meanNanos() / 1.0e6,
                unit,
                timer.calls));
    }

    private abstract static class BaseEvaluator implements GradientWrtParameterProvider {
        final String label;
        final int dim;
        final int num2x2Blocks;
        final ParameterBundle parameters;
        final CompoundParameter hmcParameter;
        final Parameter mask;
        final ProfileStats stats = new ProfileStats();
        final double[] gradient;
        final Likelihood likelihood;

        BaseEvaluator(final String label, final FixtureData data) {
            this.label = label;
            this.dim = data.dim;
            this.num2x2Blocks = data.num2x2Blocks;
            this.parameters = new ParameterBundle(label, data);
            this.hmcParameter = createHmcParameter(label, parameters);
            this.mask = createMask(parameters);
            this.gradient = new double[hmcParameter.getDimension()];
            this.likelihood = new Likelihood.Abstract(null) {
                {
                    setId(label + ".likelihood");
                }

                @Override
                protected boolean getLikelihoodKnown() {
                    return false;
                }

                @Override
                protected double calculateLogLikelihood() {
                    final long start = System.nanoTime();
                    final double logL = BaseEvaluator.this.computeLogLikelihoodRaw();
                    stats.logLikelihood.add(System.nanoTime() - start);
                    return logL;
                }
            };
        }

        @Override
        public Likelihood getLikelihood() {
            return likelihood;
        }

        @Override
        public Parameter getParameter() {
            return hmcParameter;
        }

        @Override
        public int getDimension() {
            return hmcParameter.getDimension();
        }

        @Override
        public double[] getGradientLogDensity() {
            final long start = System.nanoTime();
            Arrays.fill(gradient, 0.0);
            computeGradientRaw(gradient);
            stats.gradient.add(System.nanoTime() - start);
            return gradient;
        }

        void resetToBase() {
            parameters.resetToBase();
        }

        abstract double computeLogLikelihoodRaw();

        abstract void computeGradientRaw(double[] gradientOut);
    }

    private static final class DenseProfileEvaluator extends BaseEvaluator {
        private final DefaultTreeModel treeModel;
        private final HomogeneousCanonicalOUBranchTransitionProvider provider;
        private final SequentialCanonicalOUMessagePasser passer;
        private final CanonicalRootPrior rootPrior;
        private final double[] gradA;
        private final double[] gradQ;
        private final double[] gradMu;

        DenseProfileEvaluator(final FixtureData data) {
            super("dense homogeneous edge-rotation path", data);
            this.treeModel = buildTreeModel(data.tree);
            this.provider = new HomogeneousCanonicalOUBranchTransitionProvider(
                    treeModel,
                    new MultivariateElasticModel(parameters.blockParameter),
                    parameters.diffusionModel,
                    parameters.muParam,
                    new DefaultBranchRateModel());
            this.passer = new SequentialCanonicalOUMessagePasser(treeModel, dim);
            for (int t = 0; t < TIP_COUNT; ++t) {
                final CanonicalTipObservation observation = new CanonicalTipObservation(dim);
                observation.setObserved(data.tipValues[t]);
                passer.setTipObservation(t, observation);
            }
            this.rootPrior = new CanonicalConjugateRootPriorAdapter(
                    new ConjugateRootTraitPrior(
                            new Parameter.Default(data.beta.clone()),
                            new Parameter.Default(KAPPA0)),
                    dim);
            this.gradA = new double[parameters.nativeBlockDim + dim * dim];
            this.gradQ = new double[dim * dim];
            this.gradMu = new double[dim];
        }

        @Override
        double computeLogLikelihoodRaw() {
            return passer.computePostOrderLogLikelihood(provider, rootPrior);
        }

        @Override
        void computeGradientRaw(final double[] gradientOut) {
            long start = System.nanoTime();
            passer.computePostOrderLogLikelihood(provider, rootPrior);
            stats.phase("postorder likelihood").add(System.nanoTime() - start);

            start = System.nanoTime();
            passer.computePreOrder(provider, rootPrior);
            stats.phase("preorder messages").add(System.nanoTime() - start);

            start = System.nanoTime();
            passer.computeJointGradients(provider, gradA, gradQ, gradMu);
            stats.phase("joint gradient pullback").add(System.nanoTime() - start);

            start = System.nanoTime();
            System.arraycopy(gradA, 0, gradientOut, 0, gradA.length);
            packRowMajorAsColumnMajor(gradQ, dim, gradientOut, parameters.qOffset);
            stats.phase("gradient packing").add(System.nanoTime() - start);
        }
    }

    private static final class BlockBasisProfileEvaluator extends BaseEvaluator {
        private final ProfileTreeNode root;
        private final Scratch scratch;

        BlockBasisProfileEvaluator(final FixtureData data) {
            super("new rotated block-basis path", data);
            this.root = data.tree.deepCopy();
            allocateNodeStates(root, dim, true);
            this.scratch = new Scratch(dim, parameters, data);
        }

        @Override
        double computeLogLikelihoodRaw() {
            prepareSnapshotAndTips();
            computePostOrder(root, scratch);
            fillPriorState(scratch);
            CanonicalGaussianMessageOps.combineStates(root.selfMessage, scratch.priorState, root.nodeFullPosterior);
            return CanonicalGaussianMessageOps.normalizationShift(root.nodeFullPosterior, scratch.ws);
        }

        @Override
        void computeGradientRaw(final double[] gradientOut) {
            long start = System.nanoTime();
            prepareSnapshotAndTips();
            stats.phase("snapshot + tip rotation").add(System.nanoTime() - start);

            start = System.nanoTime();
            Arrays.fill(scratch.compressedDAcc, 0.0);
            CommonOps.fill(scratch.gSAcc, 0.0);
            CommonOps.fill(scratch.rGradAcc, 0.0);
            Arrays.fill(scratch.qGradFlatAcc, 0.0);
            computePostOrder(root, scratch);
            stats.phase("postorder in block basis").add(System.nanoTime() - start);

            start = System.nanoTime();
            fillPriorState(scratch);
            CanonicalGaussianMessageOps.combineStates(root.selfMessage, scratch.priorState, root.nodeFullPosterior);
            final double logL = CanonicalGaussianMessageOps.normalizationShift(root.nodeFullPosterior, scratch.ws);
            if (!Double.isFinite(logL)) {
                throw new AssertionError("Non-finite block-basis logL: " + logL);
            }
            stats.phase("root prior combine").add(System.nanoTime() - start);

            start = System.nanoTime();
            computeBackwardAndPullback(root, scratch.priorState, scratch);
            accumulateRootPriorVtildePullback(root.nodeFullPosterior, scratch.priorState, scratch);
            stats.phase("preorder + local pullbacks").add(System.nanoTime() - start);

            start = System.nanoTime();
            scratch.adjointBuilder.finalizeGlobalPullback(
                    scratch.snapshot,
                    scratch.gSAcc,
                    scratch.compressedDAcc,
                    scratch.rGradAcc,
                    scratch.qGradFlatAcc);
            scratch.adjointBuilder.addJacobianTermGradient(
                    scratch.snapshot,
                    TIP_COUNT,
                    scratch.rGradAcc);
            parameters.blockParameter.chainGradient(scratch.compressedDAcc, scratch.nativeGradAcc);
            stats.phase("global finalize + chain").add(System.nanoTime() - start);

            start = System.nanoTime();
            System.arraycopy(scratch.nativeGradAcc, 0, gradientOut, 0, parameters.nativeBlockDim);
            packDenseMatrixAsColumnMajor(scratch.rGradAcc, gradientOut, parameters.rotationOffset);
            packRowMajorAsColumnMajor(scratch.qGradFlatAcc, dim, gradientOut, parameters.qOffset);
            stats.phase("gradient packing").add(System.nanoTime() - start);
        }

        private void prepareSnapshotAndTips() {
            scratch.snapshot.ensureCurrent();
            for (int t = 0; t < TIP_COUNT; ++t) {
                transformToTilde(
                        scratch.snapshot,
                        scratch.tipValuesNatural[t],
                        scratch.beta,
                        dim,
                        scratch.tipValuesTilde[t]);
            }
        }
    }

    private static CompoundParameter createHmcParameter(final String label, final ParameterBundle p) {
        final CompoundParameter compound = new CompoundParameter(label + ".hmc.parameter");
        if (p.scalarParam.getDimension() > 0) {
            compound.addParameter(p.scalarParam);
        }
        compound.addParameter(p.rhoParam);
        compound.addParameter(p.thetaParam);
        compound.addParameter(p.tParam);
        compound.addParameter(p.rotationParam);
        compound.addParameter(p.qMatrix);
        return compound;
    }

    private static Parameter createMask(final ParameterBundle p) {
        final double[] mask = new double[p.hmcDimension];
        Arrays.fill(mask, 0, p.rotationOffset, 1.0);
        for (int col = 0; col < p.dim; ++col) {
            for (int row = 0; row < p.dim; ++row) {
                final int index = p.rotationOffset + col * p.dim + row;
                mask[index] = row == col ? 0.0 : 1.0;
            }
        }
        return new Parameter.Default(mask);
    }

    private static void computePostOrder(final ProfileTreeNode node, final Scratch s) {
        if (node.leaf) {
            s.builder.fillCanonicalTransition(s.snapshot, node.branchLength, node.transition);
            CanonicalGaussianMessageOps.conditionOnObservedSecondBlock(
                    node.transition,
                    s.tipValuesTilde[node.tipIndex],
                    node.pushedUpMessage);
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

    private static void computeBackwardAndPullback(final ProfileTreeNode node,
                                                   final CanonicalGaussianState outsideForNode,
                                                   final Scratch s) {
        final ProfileTreeNode left = node.left;
        final ProfileTreeNode right = node.right;

        CanonicalGaussianMessageOps.combineStates(outsideForNode, right.pushedUpMessage, s.excludingRight);
        CanonicalGaussianMessageOps.combineStates(outsideForNode, left.pushedUpMessage, s.excludingLeft);

        if (left.leaf) {
            accumulateTipPullback(node.nodeFullPosterior, s.tipValuesTilde[left.tipIndex],
                    left.transition, left.branchLength, s);
        } else {
            CanonicalGaussianMessageOps.buildPairPosterior(s.excludingRight, left.transition, left.selfMessage, s.pairState);
            CanonicalBranchMessageContributionUtils.fillFromPairState(s.pairState, s.pairWs, s.contribution);
            CanonicalTransitionAdjointUtils.fillFromCanonicalTransition(
                    left.transition, s.contribution, s.transitionAdjointWs, s.adjoints);
            s.adjointBuilder.accumulateBranchPullback(
                    s.snapshot, left.branchLength, s.adjoints, s.compressedDAcc, s.gSAcc);
            CanonicalGaussianMessageOps.pushForward(s.excludingRight, left.transition, s.ws, left.outsideMessage);
            CanonicalGaussianMessageOps.combineStates(left.selfMessage, left.outsideMessage, left.nodeFullPosterior);
        }

        if (right.leaf) {
            accumulateTipPullback(node.nodeFullPosterior, s.tipValuesTilde[right.tipIndex],
                    right.transition, right.branchLength, s);
        } else {
            CanonicalGaussianMessageOps.buildPairPosterior(s.excludingLeft, right.transition, right.selfMessage, s.pairState);
            CanonicalBranchMessageContributionUtils.fillFromPairState(s.pairState, s.pairWs, s.contribution);
            CanonicalTransitionAdjointUtils.fillFromCanonicalTransition(
                    right.transition, s.contribution, s.transitionAdjointWs, s.adjoints);
            s.adjointBuilder.accumulateBranchPullback(
                    s.snapshot, right.branchLength, s.adjoints, s.compressedDAcc, s.gSAcc);
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

    private static void accumulateTipPullback(final CanonicalGaussianState parentPosterior,
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

        CanonicalTransitionAdjointUtils.fillFromCanonicalTransition(
                transition, s.contribution, s.transitionAdjointWs, s.adjoints);
        s.adjointBuilder.accumulateBranchPullback(
                s.snapshot, dt, s.adjoints, s.compressedDAcc, s.gSAcc);
    }

    private static void accumulateRootPriorVtildePullback(final CanonicalGaussianState rootPosterior,
                                                          final CanonicalGaussianState priorState,
                                                          final Scratch s) {
        final int dim = s.dim;
        GaussianFormConverter.fillMomentsFromState(rootPosterior, s.parentMean, s.parentCovFlat, dim, s.gfcWs);

        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                s.dLogLdPriorCov[i * dim + j] =
                        -0.5 * (s.parentCovFlat[i * dim + j] + s.parentMean[i] * s.parentMean[j]);
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

    private static void fillPriorState(final Scratch s) {
        final DenseMatrix64F vTilde = s.snapshot.getStationaryCovarianceDBasis();
        for (int i = 0; i < s.dim * s.dim; ++i) {
            s.priorCovFlat[i] = vTilde.data[i] / KAPPA0;
        }
        GaussianFormConverter.fillStateFromMoments(s.zeroMean, s.priorCovFlat, s.dim, s.gfcWs, s.priorState);
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

    private static void allocateNodeStates(final ProfileTreeNode node, final int dim, final boolean isRoot) {
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

    private static void packDenseMatrixAsColumnMajor(final DenseMatrix64F matrix,
                                                     final double[] out,
                                                     final int offset) {
        final int dim = matrix.numRows;
        for (int col = 0; col < dim; ++col) {
            for (int row = 0; row < dim; ++row) {
                out[offset + col * dim + row] = matrix.get(row, col);
            }
        }
    }

    private static void packRowMajorAsColumnMajor(final double[] rowMajor,
                                                  final int dim,
                                                  final double[] out,
                                                  final int offset) {
        for (int col = 0; col < dim; ++col) {
            for (int row = 0; row < dim; ++row) {
                out[offset + col * dim + row] = rowMajor[row * dim + col];
            }
        }
    }

    private static DefaultTreeModel buildTreeModel(final ProfileTreeNode root) {
        try {
            final Tree tree = new NewickImporter(toNewick(root, true)).importTree(null);
            return new DefaultTreeModel("rotatedOUProfileTree", tree);
        } catch (IOException | Importer.ImportException e) {
            throw new RuntimeException("Failed to build profile tree model", e);
        }
    }

    private static String toNewick(final ProfileTreeNode node, final boolean isRoot) {
        final String body;
        if (node.leaf) {
            body = "taxon" + node.tipIndex;
        } else {
            body = "(" + toNewick(node.left, false) + "," + toNewick(node.right, false) + ")";
        }
        if (isRoot) {
            return body + ";";
        }
        return body + ":" + node.branchLength;
    }

    private static FixtureData buildFixtureData(final int dim) {
        final Random random = new Random(FIXTURE_SEED + dim);
        final boolean oddDim = (dim & 1) == 1;
        final int num2x2Blocks = dim / 2;

        final double[] baseRotation = new double[dim * dim];
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                baseRotation[i * dim + j] = (i == j) ? 1.0 : (random.nextDouble() - 0.5) * (0.3 / dim);
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
        final double[] baseScalar = oddDim ? new double[]{0.5 + random.nextDouble()} : new double[0];
        final double[] baseQ = randomSpdFlat(dim, random);

        final double[] beta = new double[dim];
        for (int i = 0; i < dim; ++i) {
            beta[i] = random.nextGaussian();
        }

        final double[][] tipValues = new double[TIP_COUNT][dim];
        for (int t = 0; t < TIP_COUNT; ++t) {
            for (int i = 0; i < dim; ++i) {
                tipValues[t][i] = 3.0 * Math.sin(0.13 * (t + 1) * (i + 1));
            }
        }

        final ProfileTreeNode tree = buildProfileTree(0, TIP_COUNT, random, true);
        return new FixtureData(dim, num2x2Blocks, baseRotation, baseScalar, baseRho,
                baseTheta, baseT, baseQ, beta, tipValues, tree);
    }

    private static ProfileTreeNode buildProfileTree(final int start,
                                                    final int end,
                                                    final Random random,
                                                    final boolean isRoot) {
        final ProfileTreeNode node = new ProfileTreeNode();
        if (end - start == 1) {
            node.leaf = true;
            node.tipIndex = start;
        } else {
            node.leaf = false;
            final int mid = start + (end - start) / 2;
            node.left = buildProfileTree(start, mid, random, false);
            node.right = buildProfileTree(mid, end, random, false);
        }
        if (!isRoot) {
            node.branchLength = 0.05 + random.nextDouble() * 0.55;
        }
        return node;
    }

    private static double[] randomSpdFlat(final int dim, final Random random) {
        final double[][] a = new double[dim][dim];
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                a[i][j] = random.nextDouble() - 0.5;
            }
        }
        final double[] q = new double[dim * dim];
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                double sum = 0.0;
                for (int k = 0; k < dim; ++k) {
                    sum += a[i][k] * a[j][k];
                }
                if (i == j) {
                    sum += 0.5;
                }
                q[i * dim + j] = sum;
            }
        }
        return q;
    }

    private static void fillMatrixParameter(final MatrixParameter matrix,
                                            final double[] rowMajor,
                                            final int dim) {
        for (int row = 0; row < dim; ++row) {
            for (int col = 0; col < dim; ++col) {
                matrix.setParameterValue(row, col, rowMajor[row * dim + col]);
            }
        }
    }

    private static void syncPrecisionFromVariance(final MatrixParameter qMatrix,
                                                  final MatrixParameter precisionMatrix) {
        final int dim = qMatrix.getRowDimension();
        final DenseMatrix64F variance = new DenseMatrix64F(dim, dim);
        final DenseMatrix64F precision = new DenseMatrix64F(dim, dim);
        for (int row = 0; row < dim; ++row) {
            for (int col = 0; col < dim; ++col) {
                variance.set(row, col, qMatrix.getParameterValue(row, col));
            }
        }
        if (!CommonOps.invert(variance, precision)) {
            throw new IllegalArgumentException("Variance matrix is singular in rotated OU HMC profile.");
        }
        for (int row = 0; row < dim; ++row) {
            for (int col = 0; col < dim; ++col) {
                precisionMatrix.setParameterValueQuietly(row, col, precision.get(row, col));
            }
        }
        precisionMatrix.fireParameterChangedEvent();
    }

    private static final class FixtureData {
        final int dim;
        final int num2x2Blocks;
        final double[] baseRotation;
        final double[] baseScalar;
        final double[] baseRho;
        final double[] baseTheta;
        final double[] baseT;
        final double[] baseQ;
        final double[] beta;
        final double[][] tipValues;
        final ProfileTreeNode tree;

        FixtureData(final int dim,
                    final int num2x2Blocks,
                    final double[] baseRotation,
                    final double[] baseScalar,
                    final double[] baseRho,
                    final double[] baseTheta,
                    final double[] baseT,
                    final double[] baseQ,
                    final double[] beta,
                    final double[][] tipValues,
                    final ProfileTreeNode tree) {
            this.dim = dim;
            this.num2x2Blocks = num2x2Blocks;
            this.baseRotation = baseRotation;
            this.baseScalar = baseScalar;
            this.baseRho = baseRho;
            this.baseTheta = baseTheta;
            this.baseT = baseT;
            this.baseQ = baseQ;
            this.beta = beta;
            this.tipValues = tipValues;
            this.tree = tree;
        }
    }

    private static final class ParameterBundle {
        final String label;
        final int dim;
        final int nativeBlockDim;
        final int rotationOffset;
        final int qOffset;
        final int hmcDimension;
        final double[] baseRotation;
        final double[] baseScalar;
        final double[] baseRho;
        final double[] baseTheta;
        final double[] baseT;
        final double[] baseQ;

        final MatrixParameter rotationParam;
        final Parameter scalarParam;
        final Parameter rhoParam;
        final Parameter thetaParam;
        final Parameter tParam;
        final BlockDiagonalPolarStableMatrixParameter blockParameter;
        final MatrixParameter qMatrix;
        final MatrixParameter precisionMatrix;
        final MultivariateDiffusionModel diffusionModel;
        final Parameter muParam;

        ParameterBundle(final String label, final FixtureData data) {
            this.label = label;
            this.dim = data.dim;
            this.baseRotation = data.baseRotation.clone();
            this.baseScalar = data.baseScalar.clone();
            this.baseRho = data.baseRho.clone();
            this.baseTheta = data.baseTheta.clone();
            this.baseT = data.baseT.clone();
            this.baseQ = data.baseQ.clone();

            this.rotationParam = new MatrixParameter(label + ".R", dim, dim);
            this.scalarParam = new Parameter.Default(baseScalar.length);
            this.rhoParam = new Parameter.Default(baseRho.clone());
            this.thetaParam = new Parameter.Default(baseTheta.clone());
            this.tParam = new Parameter.Default(baseT.clone());
            this.qMatrix = new MatrixParameter(label + ".Q", dim, dim);
            this.precisionMatrix = new MatrixParameter(label + ".precision", dim, dim);
            this.muParam = new Parameter.Default(data.beta.clone());

            this.blockParameter = new BlockDiagonalPolarStableMatrixParameter(
                    label + ".block",
                    rotationParam,
                    scalarParam,
                    rhoParam,
                    thetaParam,
                    tParam);
            this.diffusionModel = new MultivariateDiffusionModel(precisionMatrix);

            this.nativeBlockDim = blockParameter.getBlockDiagonalNParameters();
            this.rotationOffset = nativeBlockDim;
            this.qOffset = rotationOffset + dim * dim;
            this.hmcDimension = qOffset + dim * dim;
            resetToBase();
        }

        void resetToBase() {
            fillMatrixParameter(rotationParam, baseRotation, dim);
            fillMatrixParameter(qMatrix, baseQ, dim);
            syncPrecisionFromVariance(qMatrix, precisionMatrix);
            for (int i = 0; i < baseScalar.length; ++i) {
                scalarParam.setParameterValue(i, baseScalar[i]);
            }
            for (int i = 0; i < baseRho.length; ++i) {
                rhoParam.setParameterValue(i, baseRho[i]);
                thetaParam.setParameterValue(i, baseTheta[i]);
                tParam.setParameterValue(i, baseT[i]);
            }
        }
    }

    private static final class ProfileTreeNode {
        ProfileTreeNode left;
        ProfileTreeNode right;
        boolean leaf;
        int tipIndex;
        double branchLength;

        CanonicalGaussianTransition transition;
        CanonicalGaussianState pushedUpMessage;
        CanonicalGaussianState selfMessage;
        CanonicalGaussianState outsideMessage;
        CanonicalGaussianState nodeFullPosterior;

        ProfileTreeNode deepCopy() {
            final ProfileTreeNode copy = new ProfileTreeNode();
            copy.leaf = leaf;
            copy.tipIndex = tipIndex;
            copy.branchLength = branchLength;
            if (left != null) {
                copy.left = left.deepCopy();
            }
            if (right != null) {
                copy.right = right.deepCopy();
            }
            return copy;
        }
    }

    private static final class Scratch {
        final int dim;
        final BlockBasisGlobalSnapshot snapshot;
        final BlockBasisCanonicalTransitionBuilder builder;
        final BlockBasisCanonicalAdjointBuilder adjointBuilder;
        final double[][] tipValuesNatural;
        final double[] beta;

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

        Scratch(final int dim, final ParameterBundle parameters, final FixtureData data) {
            this.dim = dim;
            this.snapshot = new BlockBasisGlobalSnapshot(parameters.blockParameter, parameters.qMatrix);
            this.snapshot.ensureCurrent();
            this.builder = new BlockBasisCanonicalTransitionBuilder(snapshot.getStructure());
            this.adjointBuilder = new BlockBasisCanonicalAdjointBuilder(snapshot.getStructure());
            this.tipValuesNatural = data.tipValues;
            this.beta = data.beta;

            this.tipValuesTilde = new double[TIP_COUNT][dim];
            this.zeroMean = new double[dim];
            this.priorCovFlat = new double[dim * dim];
            this.priorState = new CanonicalGaussianState(dim);

            this.compressedDAcc = new double[parameters.blockParameter.getCompressedDDimension()];
            this.gSAcc = new DenseMatrix64F(dim, dim);
            this.rGradAcc = new DenseMatrix64F(dim, dim);
            this.qGradFlatAcc = new double[dim * dim];
            this.nativeGradAcc = new double[parameters.nativeBlockDim];

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

    private static final class ProfileStats {
        final Timer hmcProposal = new Timer();
        final Timer logLikelihood = new Timer();
        final Timer gradient = new Timer();
        final LinkedHashMap<String, Timer> phases = new LinkedHashMap<>();

        Timer phase(final String name) {
            Timer timer = phases.get(name);
            if (timer == null) {
                timer = new Timer();
                phases.put(name, timer);
            }
            return timer;
        }

        void reset() {
            hmcProposal.reset();
            logLikelihood.reset();
            gradient.reset();
            phases.clear();
        }
    }

    private static final class Timer {
        long nanos;
        long calls;

        void add(final long elapsedNanos) {
            nanos += elapsedNanos;
            calls++;
        }

        double meanNanos() {
            return calls == 0L ? 0.0 : nanos / (double) calls;
        }

        void reset() {
            nanos = 0L;
            calls = 0L;
        }
    }
}

package test.dr.evomodel.treedatalikelihood.continuous;

import dr.evomodel.tree.TreeModel;
import dr.evomodel.treedatalikelihood.DataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.BranchSpecificGradient;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.continuous.ContinuousTraitGradientForBranch;
import dr.evomodel.treedatalikelihood.continuous.DiffusionProcessDelegate;
import dr.evomodel.treedatalikelihood.continuous.OUDiffusionModelDelegate;
import dr.evomodel.treedatalikelihood.continuous.cdi.ContinuousDiffusionIntegrator;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.OrthogonalBlockDiagonalPolarStableMatrixParameter;
import dr.inference.model.Parameter;
import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

import java.lang.reflect.Constructor;
import java.lang.reflect.Field;
import java.lang.reflect.Method;
import java.util.List;

public class TreeCanonicalOrthogonalBlockBenchmarkTest extends TestCase {

    private static final int WARMUP_ITERS = 100;
    private static final int MEASURE_ITERS = 400;
    private static final int SCALED_BRANCH_FACTOR = 64;

    public TreeCanonicalOrthogonalBlockBenchmarkTest(final String name) {
        super(name);
    }

    public void testBenchmarkReportIsFinite() throws Exception {
        final String report = buildReport();
        System.out.println(report);
        assertTrue(report.contains("Case full-observed"));
        assertTrue(report.contains("Case partial-missing"));
        assertTrue(report.contains("Case alternate-block"));
        assertFalse(report.contains("NaN"));
        assertFalse(report.contains("Infinity"));
    }

    public static Test suite() {
        return new TestSuite(TreeCanonicalOrthogonalBlockBenchmarkTest.class);
    }

    public static void main(final String[] args) {
        junit.textui.TestRunner.run(suite());
    }

    private static String buildReport() throws Exception {
        final StringBuilder sb = new StringBuilder();
        appendCaseReport(sb, buildCase("full-observed", "buildOrthogonalBlockTreeSetups"));
        appendCaseReport(sb, buildCase("partial-missing", "buildOrthogonalBlockTreeSetupsWithPartialTipMissing"));
        appendCaseReport(sb, buildCase("alternate-block", "buildOrthogonalBlockTreeSetups", new Class<?>[]{
                double.class, double.class, double.class, double.class, double[].class
        }, new Object[]{1.1, 0.75, 0.4, 0.02, new double[]{0.35, -0.05, 0.2}}));
        return sb.toString();
    }

    private static void appendCaseReport(final StringBuilder sb, final BenchmarkCase benchmarkCase) {
        final double denseTraversalNs = timeTraversal(benchmarkCase.denseTraversalLikelihood);
        final double blockTraversalNs = timeTraversal(benchmarkCase.blockTraversalLikelihood);
        final double denseUpdateNs = timeBranchUpdate(benchmarkCase.denseTraversalLikelihood, benchmarkCase.denseDelegate, 1);
        final double blockUpdateNs = timeBranchUpdate(benchmarkCase.blockTraversalLikelihood, benchmarkCase.blockDelegate, 1);
        final double denseScaledUpdateNs = timeBranchUpdate(benchmarkCase.denseTraversalLikelihood, benchmarkCase.denseDelegate, SCALED_BRANCH_FACTOR);
        final double blockScaledUpdateNs = timeBranchUpdate(benchmarkCase.blockTraversalLikelihood, benchmarkCase.blockDelegate, SCALED_BRANCH_FACTOR);
        final double denseLikelihoodNs = timeLikelihood(benchmarkCase.denseLikelihood);
        final double blockLikelihoodNs = timeLikelihood(benchmarkCase.blockLikelihood);
        final double denseGradientNs = timeGradient(benchmarkCase.denseGradient);
        final double blockGradientNs = timeGradient(benchmarkCase.blockGradient);

        sb.append("Case ").append(benchmarkCase.label).append('\n');
        appendPair(sb, "traversal dispatch", denseTraversalNs, blockTraversalNs);
        appendPair(sb, "branch diffusion update", denseUpdateNs, blockUpdateNs);
        appendPair(sb, "branch diffusion update (x" + SCALED_BRANCH_FACTOR + " branches)", denseScaledUpdateNs, blockScaledUpdateNs);
        appendPair(sb, "tree log-likelihood", denseLikelihoodNs, blockLikelihoodNs);
        appendPair(sb, "tree analytical selection gradient", denseGradientNs, blockGradientNs);
        sb.append('\n');
    }

    private static void appendPair(final StringBuilder sb,
                                   final String label,
                                   final double denseNs,
                                   final double blockNs) {
        sb.append("  ").append(label).append('\n');
        sb.append(String.format("    dense:      %.1f ns%n", denseNs));
        sb.append(String.format("    orth-block: %.1f ns%n", blockNs));
        sb.append(String.format("    speedup:    %.2fx%n", denseNs / blockNs));
    }

    private static double timeLikelihood(final TreeDataLikelihood likelihood) {
        for (int i = 0; i < WARMUP_ITERS; ++i) {
            likelihood.makeDirty();
            likelihood.getLogLikelihood();
        }
        final long start = System.nanoTime();
        for (int i = 0; i < MEASURE_ITERS; ++i) {
            likelihood.makeDirty();
            likelihood.getLogLikelihood();
        }
        return (System.nanoTime() - start) / (double) MEASURE_ITERS;
    }

    private static double timeTraversal(final TreeDataLikelihood likelihood) {
        for (int i = 0; i < WARMUP_ITERS; ++i) {
            dispatchTraversal(likelihood);
        }
        final long start = System.nanoTime();
        for (int i = 0; i < MEASURE_ITERS; ++i) {
            dispatchTraversal(likelihood);
        }
        return (System.nanoTime() - start) / (double) MEASURE_ITERS;
    }

    private static double timeBranchUpdate(final TreeDataLikelihood likelihood,
                                           final ContinuousDataLikelihoodDelegate delegate,
                                           final int repeatFactor) {
        likelihood.makeDirty();
        likelihood.getLogLikelihood();
        for (int i = 0; i < WARMUP_ITERS; ++i) {
            runBranchUpdateOnce(likelihood, delegate, repeatFactor);
        }
        final long start = System.nanoTime();
        for (int i = 0; i < MEASURE_ITERS; ++i) {
            runBranchUpdateOnce(likelihood, delegate, repeatFactor);
        }
        return (System.nanoTime() - start) / (double) MEASURE_ITERS;
    }

    private static double timeGradient(final BranchSpecificGradient gradient) {
        for (int i = 0; i < WARMUP_ITERS; ++i) {
            gradient.getGradientLogDensity();
        }
        final long start = System.nanoTime();
        for (int i = 0; i < MEASURE_ITERS; ++i) {
            gradient.getGradientLogDensity();
        }
        return (System.nanoTime() - start) / (double) MEASURE_ITERS;
    }

    private static void runBranchUpdateOnce(final TreeDataLikelihood likelihood,
                                            final ContinuousDataLikelihoodDelegate delegate,
                                            final int repeatFactor) {
        dispatchTraversal(likelihood);
        final List<?> branchOperations = getBranchOperations(likelihood);
        final int branchCount = branchOperations.size();
        final int[] branchIndices = new int[branchCount * repeatFactor];
        final double[] branchLengths = new double[branchCount * repeatFactor];
        final double normalization = getBranchNormalization(delegate);

        for (int r = 0; r < repeatFactor; ++r) {
            final int offset = r * branchCount;
            for (int i = 0; i < branchCount; ++i) {
                final DataLikelihoodDelegate.BranchOperation op =
                        (DataLikelihoodDelegate.BranchOperation) branchOperations.get(i);
                branchIndices[offset + i] = op.getBranchNumber();
                branchLengths[offset + i] = op.getBranchLength() * normalization;
            }
        }

        final DiffusionProcessDelegate processDelegate = delegate.getDiffusionProcessDelegate();
        final ContinuousDiffusionIntegrator integrator = delegate.getIntegrator();
        final boolean flip = getFlip(delegate);
        processDelegate.updateDiffusionMatrices(integrator, branchIndices, branchLengths, branchIndices.length, flip);
    }

    private static void dispatchTraversal(final TreeDataLikelihood likelihood) {
        try {
            likelihood.makeDirty();
            final Object traversalDelegate = getFieldValue(likelihood, "treeTraversalDelegate");
            final Method dispatch = traversalDelegate.getClass().getDeclaredMethod("dispatchTreeTraversalCollectBranchAndNodeOperations");
            dispatch.setAccessible(true);
            dispatch.invoke(traversalDelegate);
        } catch (final Exception e) {
            throw new RuntimeException(e);
        }
    }

    private static BenchmarkCase buildCase(final String label, final String builderMethodName) throws Exception {
        return buildCase(label, builderMethodName, new Class<?>[0], new Object[0]);
    }

    private static BenchmarkCase buildCase(final String label,
                                           final String builderMethodName,
                                           final Class<?>[] parameterTypes,
                                           final Object[] arguments) throws Exception {
        final Class<?> testClass = Class.forName("test.dr.evomodel.treedatalikelihood.continuous.OUDiffusionKernelBridgeValidationTest");
        final Constructor<?> constructor = testClass.getConstructor(String.class);
        final Object testInstance = constructor.newInstance("benchmark");
        invokeLifecycleIfPresent(testInstance, "setUp");

        final Method builder = testClass.getDeclaredMethod(builderMethodName, parameterTypes);
        builder.setAccessible(true);
        final Object setup = builder.invoke(testInstance, arguments);

        final Object blockSetup = getFieldValue(setup, "block");
        final Object denseSetup = getFieldValue(setup, "dense");

        final TreeDataLikelihood blockLikelihood = (TreeDataLikelihood) getFieldValue(blockSetup, "treeDataLikelihood");
        final TreeDataLikelihood denseLikelihood = (TreeDataLikelihood) getFieldValue(denseSetup, "treeDataLikelihood");
        final ContinuousDataLikelihoodDelegate blockDelegate = (ContinuousDataLikelihoodDelegate) getFieldValue(blockSetup, "likelihoodDelegate");
        final ContinuousDataLikelihoodDelegate denseDelegate = (ContinuousDataLikelihoodDelegate) getFieldValue(denseSetup, "likelihoodDelegate");
        final OrthogonalBlockDiagonalPolarStableMatrixParameter blockParameter =
                (OrthogonalBlockDiagonalPolarStableMatrixParameter) getFieldValue(blockSetup, "blockParameter");
        final OUDiffusionModelDelegate denseOuDelegate =
                (OUDiffusionModelDelegate) denseDelegate.getDiffusionProcessDelegate();
        final MatrixParameterInterface denseSelectionParameter =
                denseOuDelegate.getElasticModel().getStrengthOfSelectionMatrixParameter();
        final int dimTrait = (Integer) getFieldValue(testInstance, "dimTrait");
        final TreeModel treeModel = (TreeModel) getFieldValue(testInstance, "treeModel");

        final Parameter nativeParameter = blockParameter.getParameter();
        final ContinuousTraitGradientForBranch.SelectionParameterGradient blockSelectionGradient =
                new ContinuousTraitGradientForBranch.SelectionParameterGradient(
                        dimTrait, treeModel, blockDelegate, nativeParameter, blockParameter);
        final BranchSpecificGradient blockGradient =
                new BranchSpecificGradient("trait", blockLikelihood, blockDelegate, blockSelectionGradient, nativeParameter);

        final Parameter denseParameter = flattenMatrixParameter(denseSelectionParameter);
        final ContinuousTraitGradientForBranch.SelectionParameterGradient denseSelectionGradient =
                new ContinuousTraitGradientForBranch.SelectionParameterGradient(
                        dimTrait, treeModel, denseDelegate, denseParameter, null);
        final BranchSpecificGradient denseGradient =
                new BranchSpecificGradient("trait", denseLikelihood, denseDelegate, denseSelectionGradient, denseParameter);

        return new BenchmarkCase(label, blockLikelihood, denseLikelihood, blockGradient, denseGradient,
                blockDelegate, denseDelegate, blockLikelihood, denseLikelihood);
    }

    private static Parameter flattenMatrixParameter(final MatrixParameterInterface matrixParameter) {
        final int rows = matrixParameter.getRowDimension();
        final int cols = matrixParameter.getColumnDimension();
        final Parameter.Default parameter = new Parameter.Default(rows * cols);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                parameter.setParameterValue(i * cols + j, matrixParameter.getParameterValue(i, j));
            }
        }
        return parameter;
    }

    @SuppressWarnings("unchecked")
    private static List<?> getBranchOperations(final TreeDataLikelihood likelihood) {
        try {
            final Object traversalDelegate = getFieldValue(likelihood, "treeTraversalDelegate");
            final Method getBranchOperations = traversalDelegate.getClass().getDeclaredMethod("getBranchOperations");
            getBranchOperations.setAccessible(true);
            return (List<?>) getBranchOperations.invoke(traversalDelegate);
        } catch (final Exception e) {
            throw new RuntimeException(e);
        }
    }

    private static boolean getFlip(final ContinuousDataLikelihoodDelegate delegate) {
        try {
            return (Boolean) getFieldValue(delegate, "flip");
        } catch (final Exception e) {
            throw new RuntimeException(e);
        }
    }

    private static double getBranchNormalization(final ContinuousDataLikelihoodDelegate delegate) {
        try {
            final Method getRateTransformation = delegate.getClass().getDeclaredMethod("getRateTransformation");
            getRateTransformation.setAccessible(true);
            final Object rateTransformation = getRateTransformation.invoke(delegate);
            final Method getNormalization = rateTransformation.getClass().getMethod("getNormalization");
            return (Double) getNormalization.invoke(rateTransformation);
        } catch (final Exception e) {
            throw new RuntimeException(e);
        }
    }

    private static Object getFieldValue(final Object target, final String fieldName) throws Exception {
        Class<?> type = target.getClass();
        while (type != null) {
            try {
                final Field field = type.getDeclaredField(fieldName);
                field.setAccessible(true);
                return field.get(target);
            } catch (final NoSuchFieldException ignored) {
                type = type.getSuperclass();
            }
        }
        throw new NoSuchFieldException(fieldName);
    }

    private static void invokeLifecycleIfPresent(final Object target, final String methodName) throws Exception {
        Class<?> type = target.getClass();
        while (type != null) {
            try {
                final Method method = type.getDeclaredMethod(methodName);
                method.setAccessible(true);
                method.invoke(target);
                return;
            } catch (final NoSuchMethodException ignored) {
                type = type.getSuperclass();
            }
        }
    }

    private static final class BenchmarkCase {
        final String label;
        final TreeDataLikelihood blockLikelihood;
        final TreeDataLikelihood denseLikelihood;
        final BranchSpecificGradient blockGradient;
        final BranchSpecificGradient denseGradient;
        final ContinuousDataLikelihoodDelegate blockDelegate;
        final ContinuousDataLikelihoodDelegate denseDelegate;
        final TreeDataLikelihood blockTraversalLikelihood;
        final TreeDataLikelihood denseTraversalLikelihood;

        private BenchmarkCase(final String label,
                              final TreeDataLikelihood blockLikelihood,
                              final TreeDataLikelihood denseLikelihood,
                              final BranchSpecificGradient blockGradient,
                              final BranchSpecificGradient denseGradient,
                              final ContinuousDataLikelihoodDelegate blockDelegate,
                              final ContinuousDataLikelihoodDelegate denseDelegate,
                              final TreeDataLikelihood blockTraversalLikelihood,
                              final TreeDataLikelihood denseTraversalLikelihood) {
            this.label = label;
            this.blockLikelihood = blockLikelihood;
            this.denseLikelihood = denseLikelihood;
            this.blockGradient = blockGradient;
            this.denseGradient = denseGradient;
            this.blockDelegate = blockDelegate;
            this.denseDelegate = denseDelegate;
            this.blockTraversalLikelihood = blockTraversalLikelihood;
            this.denseTraversalLikelihood = denseTraversalLikelihood;
        }
    }
}

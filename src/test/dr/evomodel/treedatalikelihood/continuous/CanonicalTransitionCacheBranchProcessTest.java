package test.dr.evomodel.treedatalikelihood.continuous;

import dr.evomodel.continuous.ou.OUProcessModel;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.CanonicalPreparedBranchSnapshot;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.message.CanonicalGaussianTransition;
import dr.inference.model.MatrixParameter;
import dr.inference.model.Parameter;
import junit.framework.TestCase;

import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationHandler;
import java.lang.reflect.Method;
import java.lang.reflect.Proxy;
import java.util.Arrays;

public final class CanonicalTransitionCacheBranchProcessTest extends TestCase {

    public void testCacheDispatchesOneSegmentBranchesByProcessIndex() throws Exception {
        final OUProcessModel[] processModels = {
                makeProcess("p0", 0.25, 0.0, 1.0),
                makeProcess("p1", 1.40, 1.5, 0.7)
        };
        final double[] branchTimes = {0.0, 0.6, 0.6};
        final int[] branchProcessIndices = {0, 0, 1};

        final Object provider = fakeBranchProcessProvider(
                processModels,
                branchTimes,
                branchProcessIndices);
        final Object cache = newTransitionCache(provider);

        final CanonicalGaussianTransition actualProcess0 = new CanonicalGaussianTransition(1);
        final CanonicalGaussianTransition actualProcess1 = new CanonicalGaussianTransition(1);
        fillTransition(cache, 1, actualProcess0);
        fillTransition(cache, 2, actualProcess1);

        final CanonicalGaussianTransition expectedProcess0 = new CanonicalGaussianTransition(1);
        final CanonicalGaussianTransition expectedProcess1 = new CanonicalGaussianTransition(1);
        processModels[0].getCanonicalKernel().fillCanonicalTransition(branchTimes[1], expectedProcess0);
        processModels[1].getCanonicalKernel().fillCanonicalTransition(branchTimes[2], expectedProcess1);

        assertTransitionEquals("child 1 must use process 0", expectedProcess0, actualProcess0);
        assertTransitionEquals("child 2 must use process 1", expectedProcess1, actualProcess1);
        assertFalse("process-specific transitions should differ",
                close(expectedProcess0.precisionYY[0], actualProcess1.precisionYY[0]));

        assertEquals(0, getSnapshot(cache, 1).getProcessIndex());
        assertEquals(1, getSnapshot(cache, 2).getProcessIndex());
    }

    private static Object fakeBranchProcessProvider(final OUProcessModel[] processModels,
                                                    final double[] branchTimes,
                                                    final int[] branchProcessIndices) throws Exception {
        final Class<?> providerInterface = Class.forName(
                "dr.evomodel.treedatalikelihood.continuousNew.canonical.adapter."
                + "CanonicalOUBranchProcessProvider");
        final InvocationHandler handler = (proxy, method, args) -> {
            final String name = method.getName();
            if ("toString".equals(name)) {
                return "fake branch-process provider";
            }
            if ("hashCode".equals(name)) {
                return System.identityHashCode(proxy);
            }
            if ("equals".equals(name)) {
                return proxy == args[0];
            }
            if ("getDimension".equals(name)) {
                return 1;
            }
            if ("getNodeCount".equals(name)) {
                return branchTimes.length;
            }
            if ("getProcessModelCount".equals(name)) {
                return processModels.length;
            }
            if ("getProcessModel".equals(name)) {
                if (args == null || args.length == 0) {
                    return processModels[0];
                }
                return processModels[(Integer) args[0]];
            }
            if ("getSegmentCount".equals(name)) {
                final int child = (Integer) args[0];
                return child == 0 ? 0 : 1;
            }
            if ("getSegmentProcessIndex".equals(name)) {
                return branchProcessIndices[(Integer) args[0]];
            }
            if ("getSegmentEffectiveTime".equals(name)
                    || "getEffectiveBranchLength".equals(name)) {
                return branchTimes[(Integer) args[0]];
            }
            if ("fillCanonicalTransitionForLength".equals(name)) {
                processModels[0].getCanonicalKernel()
                        .fillCanonicalTransition((Double) args[0], (CanonicalGaussianTransition) args[1]);
                return null;
            }
            if ("createNativeBranchMessageProvider".equals(name)) {
                return null;
            }
            if ("registerWith".equals(name)
                    || "refreshSnapshot".equals(name)
                    || "beginPhase".equals(name)
                    || "endPhase".equals(name)) {
                return null;
            }
            throw new UnsupportedOperationException("Unexpected fake provider call: " + method);
        };
        return Proxy.newProxyInstance(
                providerInterface.getClassLoader(),
                new Class<?>[]{providerInterface},
                handler);
    }

    private static Object newTransitionCache(final Object provider) throws Exception {
        final Class<?> providerInterface = Class.forName(
                "dr.evomodel.treedatalikelihood.continuousNew.canonical.adapter."
                + "CanonicalOUBranchProcessProvider");
        final Class<?> optionsClass = Class.forName(
                "dr.evomodel.treedatalikelihood.continuousNew.canonical.adapter."
                + "CanonicalTransitionCacheOptions");
        final Method optionsFactory = optionsClass.getDeclaredMethod("fromSystemProperties");
        optionsFactory.setAccessible(true);
        final Object options = optionsFactory.invoke(null);

        final Class<?> cacheClass = Class.forName(
                "dr.evomodel.treedatalikelihood.continuousNew.canonical.adapter."
                + "CanonicalTransitionCache");
        final Constructor<?> constructor =
                cacheClass.getDeclaredConstructor(providerInterface, optionsClass);
        constructor.setAccessible(true);
        return constructor.newInstance(provider, options);
    }

    private static void fillTransition(final Object cache,
                                       final int childNodeIndex,
                                       final CanonicalGaussianTransition out) throws Exception {
        final Method method = cache.getClass().getDeclaredMethod(
                "fillTransition",
                int.class,
                CanonicalGaussianTransition.class);
        method.setAccessible(true);
        method.invoke(cache, childNodeIndex, out);
    }

    private static CanonicalPreparedBranchSnapshot getSnapshot(final Object cache,
                                                              final int childNodeIndex) throws Exception {
        final Method method = cache.getClass().getDeclaredMethod(
                "getPreparedBranchSnapshot",
                int.class);
        method.setAccessible(true);
        return (CanonicalPreparedBranchSnapshot) method.invoke(cache, childNodeIndex);
    }

    private static OUProcessModel makeProcess(final String label,
                                              final double selection,
                                              final double mean,
                                              final double diffusion) {
        return new OUProcessModel(
                label,
                1,
                makeMatrixParameter(label + ".A", new double[][]{{selection}}),
                makeMatrixParameter(label + ".Q", new double[][]{{diffusion}}),
                new Parameter.Default(label + ".mu", new double[]{mean}),
                makeMatrixParameter(label + ".P0", new double[][]{{1.0}}),
                OUProcessModel.CovarianceGradientMethod.VAN_LOAN_ADJOINT);
    }

    private static MatrixParameter makeMatrixParameter(final String name,
                                                       final double[][] values) {
        final int rows = values.length;
        final int cols = values[0].length;
        final Parameter[] columns = new Parameter[cols];
        for (int j = 0; j < cols; ++j) {
            final double[] column = new double[rows];
            for (int i = 0; i < rows; ++i) {
                column[i] = values[i][j];
            }
            columns[j] = new Parameter.Default(name + ".col" + j, column);
        }
        return new MatrixParameter(name, columns);
    }

    private static void assertTransitionEquals(final String label,
                                               final CanonicalGaussianTransition expected,
                                               final CanonicalGaussianTransition actual) {
        assertArrayEquals(label + " precisionXX", expected.precisionXX, actual.precisionXX);
        assertArrayEquals(label + " precisionXY", expected.precisionXY, actual.precisionXY);
        assertArrayEquals(label + " precisionYX", expected.precisionYX, actual.precisionYX);
        assertArrayEquals(label + " precisionYY", expected.precisionYY, actual.precisionYY);
        assertArrayEquals(label + " informationX", expected.informationX, actual.informationX);
        assertArrayEquals(label + " informationY", expected.informationY, actual.informationY);
        assertEquals(label + " logNormalizer", expected.logNormalizer, actual.logNormalizer, 1.0e-12);
    }

    private static void assertArrayEquals(final String label,
                                          final double[] expected,
                                          final double[] actual) {
        assertEquals(label + " length", expected.length, actual.length);
        for (int i = 0; i < expected.length; ++i) {
            assertEquals(label + "[" + i + "] expected=" + Arrays.toString(expected)
                            + " actual=" + Arrays.toString(actual),
                    expected[i],
                    actual[i],
                    1.0e-12);
        }
    }

    private static boolean close(final double expected, final double actual) {
        return Math.abs(expected - actual) <= 1.0e-12;
    }
}

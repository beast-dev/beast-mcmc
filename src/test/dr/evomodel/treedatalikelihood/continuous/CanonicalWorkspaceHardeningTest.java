package test.dr.evomodel.treedatalikelihood.continuous;

import junit.framework.TestCase;

import java.lang.reflect.Constructor;
import java.lang.reflect.Field;
import java.lang.reflect.Method;
import java.util.EnumSet;

public final class CanonicalWorkspaceHardeningTest extends TestCase {

    public void testTransitionMatrixCacheInvalidatesAcrossPasses() throws Exception {
        final Class<?> workspaceClass =
                Class.forName("dr.evomodel.treedatalikelihood.continuous.canonical.GradientPullbackWorkspace");
        final Object workspace = newInstance(workspaceClass, int.class, 2);
        final Method hasTransitionMatrix =
                method(workspaceClass, "hasTransitionMatrix", double.class);
        final Method cacheTransitionMatrix =
                method(workspaceClass, "cacheTransitionMatrix", double.class);
        final Method invalidateTransitionMatrixCache =
                method(workspaceClass, "invalidateTransitionMatrixCache");
        final Method clearLocalGradientBuffers =
                method(workspaceClass, "clearLocalGradientBuffers", int.class, int.class, int.class);

        assertFalse((Boolean) hasTransitionMatrix.invoke(workspace, 0.5));
        cacheTransitionMatrix.invoke(workspace, 0.5);
        assertTrue((Boolean) hasTransitionMatrix.invoke(workspace, 0.5));
        assertFalse((Boolean) hasTransitionMatrix.invoke(workspace, 0.75));

        invalidateTransitionMatrixCache.invoke(workspace);
        assertFalse((Boolean) hasTransitionMatrix.invoke(workspace, 0.5));
        cacheTransitionMatrix.invoke(workspace, 0.5);
        assertTrue((Boolean) hasTransitionMatrix.invoke(workspace, 0.5));

        clearLocalGradientBuffers.invoke(workspace, 4, 2, 2);
        assertFalse((Boolean) hasTransitionMatrix.invoke(workspace, 0.5));
    }

    public void testWorkspaceFactoryRejectsInvalidDimension() throws Exception {
        final Class<?> factoryClass =
                Class.forName("dr.evomodel.treedatalikelihood.continuous.canonical.WorkspaceFactory");
        final Method factory = method(factoryClass, "branchGradientWorkspace", int.class);
        try {
            factory.invoke(null, 0);
            fail("Expected invalid dimension to be rejected");
        } catch (java.lang.reflect.InvocationTargetException expected) {
            assertTrue(expected.getCause() instanceof IllegalArgumentException);
        }
    }

    public void testWorkspaceFactoryConstructsValidatedWorkspace() throws Exception {
        final int[] dimensions = {1, 2, 3, 8, 16};
        final Class<?> factoryClass =
                Class.forName("dr.evomodel.treedatalikelihood.continuous.canonical.WorkspaceFactory");
        final Method factory = method(factoryClass, "branchGradientWorkspace", int.class);
        for (int dim : dimensions) {
            final Object workspace = factory.invoke(null, dim);
            final Class<?> workspaceClass = workspace.getClass();

            assertEquals(dim, ((double[]) field(workspaceClass, "mean").get(workspace)).length);
            assertEquals(dim * dim, ((double[]) field(workspaceClass, "transitionMatrixFlat").get(workspace)).length);
            assertEquals(dim * dim, ((double[]) field(workspaceClass, "localGradientA").get(workspace)).length);
            assertEquals(1, ((double[]) field(workspaceClass, "localGradientMuScalar").get(workspace)).length);

            method(workspaceClass, "validate", int.class).invoke(workspace, dim);
        }
    }

    public void testWorkspaceFactoryAllocatesOnlyRequestedCapabilities() throws Exception {
        final Class<?> factoryClass =
                Class.forName("dr.evomodel.treedatalikelihood.continuous.canonical.WorkspaceFactory");
        final Class<?> capabilityClass =
                Class.forName("dr.evomodel.treedatalikelihood.continuous.canonical.WorkspaceCapability");
        final Method factory = method(factoryClass, "branchGradientWorkspace", int.class, EnumSet.class);

        final Object traversalOnly = factory.invoke(
                null,
                3,
                capabilitySet(capabilityClass, "TRAVERSAL"));
        assertNotNull(field(traversalOnly.getClass(), "traversal").get(traversalOnly));
        assertNull(field(traversalOnly.getClass(), "adjoint").get(traversalOnly));
        assertNull(field(traversalOnly.getClass(), "gradient").get(traversalOnly));

        final Object denseGradientOnly = factory.invoke(
                null,
                3,
                capabilitySet(capabilityClass, "DENSE_GRADIENT"));
        assertNull(field(denseGradientOnly.getClass(), "traversal").get(denseGradientOnly));
        assertNull(field(denseGradientOnly.getClass(), "adjoint").get(denseGradientOnly));
        final Object gradient = field(denseGradientOnly.getClass(), "gradient").get(denseGradientOnly);
        assertNotNull(gradient);
        assertNotNull(field(gradient.getClass(), "dense").get(gradient));
        assertNull(field(gradient.getClass(), "orthogonal").get(gradient));

        final Object fullGradient = factory.invoke(
                null,
                3,
                capabilitySet(capabilityClass, "DENSE_GRADIENT", "ORTHOGONAL_BLOCK_GRADIENT"));
        final Object fullGradientWorkspace = field(fullGradient.getClass(), "gradient").get(fullGradient);
        assertNotNull(field(fullGradientWorkspace.getClass(), "dense").get(fullGradientWorkspace));
        assertNotNull(field(fullGradientWorkspace.getClass(), "orthogonal").get(fullGradientWorkspace));
    }

    public void testWorkspaceFactoryRejectsUnsupportedCapabilityCombinations() throws Exception {
        final Class<?> factoryClass =
                Class.forName("dr.evomodel.treedatalikelihood.continuous.canonical.WorkspaceFactory");
        final Class<?> capabilityClass =
                Class.forName("dr.evomodel.treedatalikelihood.continuous.canonical.WorkspaceCapability");
        final Method factory = method(factoryClass, "branchGradientWorkspace", int.class, EnumSet.class);

        assertIllegalArgument(factory, capabilitySet(capabilityClass, "PARTIAL_OBSERVATION"));
        assertIllegalArgument(factory, capabilitySet(capabilityClass, "ORTHOGONAL_BLOCK_GRADIENT"));
    }

    private static Object newInstance(final Class<?> clazz,
                                      final Class<?> parameterType,
                                      final Object argument) throws Exception {
        final Constructor<?> constructor = clazz.getDeclaredConstructor(parameterType);
        constructor.setAccessible(true);
        return constructor.newInstance(argument);
    }

    private static Method method(final Class<?> clazz,
                                 final String name,
                                 final Class<?>... parameterTypes) throws Exception {
        final Method method = clazz.getDeclaredMethod(name, parameterTypes);
        method.setAccessible(true);
        return method;
    }

    private static Field field(final Class<?> clazz, final String name) throws Exception {
        final Field field = clazz.getDeclaredField(name);
        field.setAccessible(true);
        return field;
    }

    @SuppressWarnings({"unchecked", "rawtypes"})
    private static EnumSet<?> capabilitySet(final Class<?> capabilityClass,
                                            final String... names) {
        final EnumSet capabilities = EnumSet.noneOf((Class<? extends Enum>) capabilityClass);
        for (String name : names) {
            capabilities.add(Enum.valueOf((Class<? extends Enum>) capabilityClass, name));
        }
        return capabilities;
    }

    private static void assertIllegalArgument(final Method factory,
                                              final EnumSet<?> capabilities) throws Exception {
        try {
            factory.invoke(null, 3, capabilities);
            fail("Expected unsupported capability combination to be rejected");
        } catch (java.lang.reflect.InvocationTargetException expected) {
            assertTrue(expected.getCause() instanceof IllegalArgumentException);
        }
    }
}

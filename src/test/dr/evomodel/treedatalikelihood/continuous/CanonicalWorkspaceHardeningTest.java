package test.dr.evomodel.treedatalikelihood.continuous;

import dr.evomodel.continuous.ou.DenseSelectionMatrixParameterization;
import dr.evomodel.continuous.ou.SelectionMatrixParameterization;
import dr.evomodel.continuous.ou.SelectionMatrixParameterizationFactory;
import dr.evomodel.continuous.ou.diagonal.DiagonalSelectionMatrixParameterization;
import dr.inference.model.DiagonalMatrix;
import dr.inference.model.Parameter;
import junit.framework.TestCase;

import java.lang.reflect.Constructor;
import java.lang.reflect.Field;
import java.lang.reflect.Method;
import java.io.File;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.EnumSet;

public final class CanonicalWorkspaceHardeningTest extends TestCase {

    public void testTransitionMatrixCacheInvalidatesAcrossPasses() throws Exception {
        final Class<?> workspaceClass =
                Class.forName("dr.evomodel.treedatalikelihood.continuousNew.canonical.workspace.GradientPullbackWorkspace");
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
                Class.forName("dr.evomodel.treedatalikelihood.continuousNew.canonical.workspace.WorkspaceFactory");
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
                Class.forName("dr.evomodel.treedatalikelihood.continuousNew.canonical.workspace.WorkspaceFactory");
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
                Class.forName("dr.evomodel.treedatalikelihood.continuousNew.canonical.workspace.WorkspaceFactory");
        final Class<?> capabilityClass =
                Class.forName("dr.evomodel.treedatalikelihood.continuousNew.canonical.workspace.WorkspaceCapability");
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
        assertNull(field(gradient.getClass(), "specialized").get(gradient));

        final Object fullGradient = factory.invoke(
                null,
                3,
                capabilitySet(capabilityClass, "DENSE_GRADIENT", "SPECIALIZED_GRADIENT"));
        final Object fullGradientWorkspace = field(fullGradient.getClass(), "gradient").get(fullGradient);
        assertNotNull(field(fullGradientWorkspace.getClass(), "dense").get(fullGradientWorkspace));
        assertNotNull(field(fullGradientWorkspace.getClass(), "specialized").get(fullGradientWorkspace));
    }

    public void testWorkspaceFactoryRejectsUnsupportedCapabilityCombinations() throws Exception {
        final Class<?> factoryClass =
                Class.forName("dr.evomodel.treedatalikelihood.continuousNew.canonical.workspace.WorkspaceFactory");
        final Class<?> capabilityClass =
                Class.forName("dr.evomodel.treedatalikelihood.continuousNew.canonical.workspace.WorkspaceCapability");
        final Method factory = method(factoryClass, "branchGradientWorkspace", int.class, EnumSet.class);

        assertIllegalArgument(factory, capabilitySet(capabilityClass, "PARTIAL_OBSERVATION"));
        assertIllegalArgument(factory, capabilitySet(capabilityClass, "SPECIALIZED_GRADIENT"));
    }

    public void testDiagonalSelectionParameterizationIsSiblingNotDenseSubclass() throws Exception {
        final File diagonalSelection = sourceFile(
                "src/dr/evomodel/continuous/ou/diagonal/DiagonalSelectionMatrixParameterization.java");
        final String source = read(diagonalSelection);
        assertFalse(source.contains("extends DenseSelectionMatrixParameterization"));

        final DiagonalMatrix diagonalMatrix =
                new DiagonalMatrix(new Parameter.Default(new double[]{0.35, 0.9, 1.7}));
        final SelectionMatrixParameterization selection =
                SelectionMatrixParameterizationFactory.create(diagonalMatrix);
        assertTrue(selection instanceof DiagonalSelectionMatrixParameterization);
        assertFalse(selection instanceof DenseSelectionMatrixParameterization);
    }

    public void testDiagonalSelectionTransitionGradientMatchesDenseDiagonalEntries() {
        final DiagonalMatrix diagonalMatrix =
                new DiagonalMatrix(new Parameter.Default(new double[]{0.35, 0.9, 1.7}));
        final DiagonalSelectionMatrixParameterization diagonalSelection =
                new DiagonalSelectionMatrixParameterization(diagonalMatrix);
        final DenseSelectionMatrixParameterization denseSelection =
                new DenseSelectionMatrixParameterization(diagonalMatrix);
        final int dim = 3;
        final double dt = 0.47;
        final double[] stationaryMean = new double[]{0.2, -0.4, 0.6};
        final double[] dLogL_dF = new double[]{
                0.11, -0.07, 0.03,
                0.05, -0.19, 0.13,
                -0.17, 0.23, 0.29
        };
        final double[] dLogL_df = new double[]{0.31, -0.41, 0.53};
        final double[] diagonalGradient = new double[dim * dim];
        final double[] denseGradient = new double[dim * dim];

        diagonalSelection.accumulateGradientFromTransitionFlat(
                dt, stationaryMean, dLogL_dF, dLogL_df, diagonalGradient);
        denseSelection.accumulateGradientFromTransitionFlat(
                dt, stationaryMean, dLogL_dF, dLogL_df, denseGradient);

        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                final int offset = i * dim + j;
                if (i == j) {
                    assertEquals(denseGradient[offset], diagonalGradient[offset], 1.0e-10);
                } else {
                    assertEquals(0.0, diagonalGradient[offset], 0.0);
                }
            }
        }
    }

    public void testDiagonalOuPackageDoesNotDependOnTreeObservationStorage() throws Exception {
        final File diagonalDir = sourceFile("src/dr/evomodel/continuous/ou/diagonal");
        assertTrue("Missing diagonal OU source directory: " + diagonalDir, diagonalDir.isDirectory());
        assertNoForbiddenText(
                diagonalDir,
                "dr.evomodel.treedatalikelihood.continuousNew.observationmodel",
                "dr.evomodel.treedatalikelihood.continuousNew.canonical.message",
                "BranchGradientWorkspace",
                "CanonicalNativeBranchMessageProvider");
    }

    public void testHomogeneousOuProviderDoesNotKnowDiagonalMatrixType() throws Exception {
        final File provider = sourceFile(
                "src/dr/evomodel/treedatalikelihood/continuousNew/canonical/adapter/"
                + "HomogeneousCanonicalOUBranchTransitionProvider.java");
        final String source = read(provider);
        assertFalse(source.contains("DiagonalMatrix"));
    }

    public void testHomogeneousOuProviderDelegatesBranchProcessDetails() throws Exception {
        final File provider = sourceFile(
                "src/dr/evomodel/treedatalikelihood/continuousNew/canonical/adapter/"
                + "HomogeneousCanonicalOUBranchTransitionProvider.java");
        final String source = read(provider);
        assertTrue(source.contains("CanonicalOUBranchProcessProvider"));
        assertFalse(source.contains("CanonicalOUDiffusionSnapshot"));
        assertFalse(source.contains("CanonicalBranchLengthSnapshot"));
        assertFalse(source.contains("CanonicalNativeBranchMessageAdapterFactory"));
        assertFalse(source.contains("new OUProcessModel"));
    }

    public void testCanonicalTransitionCacheIsBranchProcessAwareWithFastPath() throws Exception {
        final File cache = sourceFile(
                "src/dr/evomodel/treedatalikelihood/continuousNew/canonical/adapter/"
                + "CanonicalTransitionCache.java");
        final String source = read(cache);
        assertTrue(source.contains("CanonicalOUBranchProcessProvider"));
        assertTrue(source.contains("homogeneousFastPath"));
        assertTrue(source.contains("isHomogeneousSingleSegment"));
        assertTrue(source.contains("getSegmentProcessIndex"));
        assertTrue(source.contains("getSegmentEffectiveTime"));
        assertFalse(source.contains("interface BranchLengthProvider"));
    }

    public void testGradientInputsCarryBranchProcessIdentity() throws Exception {
        final String inputs = read(sourceFile(
                "src/dr/evomodel/treedatalikelihood/continuousNew/canonical/gradient/"
                + "BranchGradientInputs.java"));
        assertTrue(inputs.contains("activeProcessIndices"));
        assertTrue(inputs.contains("getProcessIndex"));

        final String preparer = read(sourceFile(
                "src/dr/evomodel/treedatalikelihood/continuousNew/canonical/gradient/"
                + "CanonicalBranchAdjointPreparer.java"));
        assertFalse(preparer.contains("getProcessModel()"));

        final String treeGradient = read(sourceFile(
                "src/dr/evomodel/treedatalikelihood/continuousNew/canonical/gradient/"
                + "CanonicalTreeGradientEngine.java"));
        assertTrue(treeGradient.contains("inputs.getProcessIndex"));
        assertTrue(treeGradient.contains("validateSharedGradientTarget"));

        final String branchLengthGradient = read(sourceFile(
                "src/dr/evomodel/treedatalikelihood/continuousNew/canonical/gradient/"
                + "CanonicalBranchLengthGradientEngine.java"));
        assertTrue(branchLengthGradient.contains("inputs.getProcessIndex"));
    }

    public void testNativeMessageAdapterFactoryUsesKernelCapabilities() throws Exception {
        final File factory = sourceFile(
                "src/dr/evomodel/treedatalikelihood/continuousNew/canonical/adapter/"
                + "CanonicalNativeBranchMessageAdapterFactory.java");
        final String source = read(factory);
        assertFalse(source.contains("dr.evomodel.continuous.ou.diagonal"));
        assertFalse(source.contains("instanceof DiagonalOUBranchKernel"));
        assertTrue(source.contains("DiagonalOUBranchTransformKernel"));
        assertTrue(source.contains("getCapability"));
    }

    public void testObservationModelsDoNotDependOnFullNativeProvider() throws Exception {
        final File observationDir = sourceFile(
                "src/dr/evomodel/treedatalikelihood/continuousNew/observationmodel");
        assertTrue("Missing observation-model source directory: " + observationDir,
                observationDir.isDirectory());
        assertNoForbiddenText(observationDir, "CanonicalNativeBranchMessageProvider");
    }

    public void testNativeTransitionPreloadPolicyIsTraversalAware() throws Exception {
        final File runner = sourceFile(
                "src/dr/evomodel/treedatalikelihood/continuousNew/canonical/CanonicalTraversalRunner.java");
        final String runnerSource = read(runner);
        assertFalse(runnerSource.contains("shouldSkipCanonicalTransitionPreload"));
        assertTrue(runnerSource.contains("postorderNativeMessagesCoverTraversal"));

        final File adapter = sourceFile(
                "src/dr/evomodel/treedatalikelihood/continuousNew/canonical/adapter/"
                + "DiagonalOUCanonicalBranchMessageAdapter.java");
        final String adapterSource = read(adapter);
        assertFalse(adapterSource.contains("shouldSkipCanonicalTransitionPreload"));
        assertTrue(adapterSource.contains("supportsPostorderTipMessage"));
    }

    public void testNativeDiagonalMessagesAreNotControlledBySystemProperty() throws Exception {
        final File continuousDir = sourceFile(
                "src/dr/evomodel/treedatalikelihood/continuous");
        assertNoForbiddenText(continuousDir, "canonicalUseNativeDiagonalOUMessages");
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

    private static File sourceFile(final String relativePath) {
        return new File(System.getProperty("user.dir"), relativePath);
    }

    private static void assertNoForbiddenText(final File file,
                                              final String... forbidden) throws Exception {
        if (file.isDirectory()) {
            final File[] children = file.listFiles();
            assertNotNull(children);
            for (File child : children) {
                if (child.isDirectory() || child.getName().endsWith(".java")) {
                    assertNoForbiddenText(child, forbidden);
                }
            }
            return;
        }
        final String source = read(file);
        for (String text : forbidden) {
            assertFalse(file + " must not contain " + text, source.contains(text));
        }
    }

    private static String read(final File file) throws Exception {
        return new String(Files.readAllBytes(file.toPath()), StandardCharsets.UTF_8);
    }
}

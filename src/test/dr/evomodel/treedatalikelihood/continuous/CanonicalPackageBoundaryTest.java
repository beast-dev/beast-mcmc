package test.dr.evomodel.treedatalikelihood.continuous;

import junit.framework.TestCase;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

public final class CanonicalPackageBoundaryTest extends TestCase {

    private static final String BLOCK_DIAGONAL_BACKEND_PACKAGE =
            "dr.evomodel.continuous.ou.blockdiagonal";
    private static final String WILDCARD_IMPORT = "wildcard import";
    private static final String SYSTEM_PRINT = "system print";

    public void testCanonicalArchitectureNoteExists() {
        final File note = new File("docs/canonical-architecture.md");
        assertTrue("Missing canonical architecture note: " + note.getPath(), note.isFile());
    }

    public void testCanonicalArchitectureNoteLocksGradientOwnership() throws IOException {
        final String note = readFile(new File("docs/canonical-architecture.md"));
        assertTrue("Architecture note must identify canonical gradient adapter ownership",
                note.contains("CanonicalOUGradientAdapter"));
        assertTrue("Architecture note must identify canonical integrator ownership",
                note.contains("CanonicalOUTreeLikelihoodIntegrator"));
        assertTrue("Architecture note must identify reusable branch gradient inputs",
                note.contains("BranchGradientInputs"));
        assertTrue("Architecture note must keep BranchSpecificGradient out of canonical backend ownership",
                note.contains("BranchSpecificGradient"));
        assertTrue("Architecture note must keep ContinuousTraitGradientForBranch out of canonical backend ownership",
                note.contains("ContinuousTraitGradientForBranch"));
    }

    public void testCanonicalArchitectureNoteLocksBlockDiagonalPullbackContract() throws IOException {
        final String note = readFile(new File("docs/canonical-architecture.md"));
        assertTrue("Architecture note must identify the specialized selection pullback path",
                note.contains("SpecializedCanonicalSelectionGradientPullback"));
        assertTrue("Architecture note must require fill/accumulate APIs",
                note.contains("fill/accumulate APIs"));
        assertTrue("Architecture note must require row-major buffers",
                note.contains("flat row-major"));
        assertTrue("Architecture note must reject double[][] in canonical gradient internals",
                note.contains("Do not introduce `double[][]` into canonical gradient internals"));
    }

    public void testCanonicalTreeLayerDoesNotImportBlockDiagonalBackend() throws IOException {
        assertNoSourceMatch(
                new File("src/dr/evomodel/treedatalikelihood/continuousNew/canonical"),
                BLOCK_DIAGONAL_BACKEND_PACKAGE);
        assertNoSourceMatch(
                new File("src/dr/evomodel/treedatalikelihood/continuous"),
                BLOCK_DIAGONAL_BACKEND_PACKAGE,
                "OUCanonical");
    }

    public void testCanonicalMachineryDoesNotUseDirectSystemPrints() throws IOException {
        assertNoSourceMatch(
                new File("src/dr/evomodel/treedatalikelihood/continuousNew/canonical"),
                SYSTEM_PRINT);
        assertNoSourceMatch(
                new File("src/dr/evomodel/continuous/canonical"),
                SYSTEM_PRINT);
        assertNoSourceMatch(
                new File("src/dr/evomodel/continuous/ou/canonical"),
                SYSTEM_PRINT);
        assertNoSourceMatch(
                new File("src/dr/evomodel/continuous/ou/blockdiagonal"),
                SYSTEM_PRINT);
        assertNoSourceMatch(
                new File("src/dr/evomodel/treedatalikelihood/continuous"),
                SYSTEM_PRINT,
                "Canonical");
    }

    public void testParentContinuousPackageDoesNotGainCanonicalImplementationClasses() {
        final File root = new File("src/dr/evomodel/treedatalikelihood/continuous");
        final File[] children = root.listFiles();
        assertNotNull("Source root does not exist: " + root.getPath(), children);
        for (File child : children) {
            if (child.isFile() && child.getName().startsWith("OUCanonical") && child.getName().endsWith(".java")) {
                fail("Parent continuous package contains canonical implementation class: " + child.getPath());
            }
        }
    }

    public void testCanonicalMachineryDoesNotUseWildcardImports() throws IOException {
        assertNoSourceMatch(
                new File("src/dr/evomodel/treedatalikelihood/continuousNew/canonical"),
                WILDCARD_IMPORT);
        assertNoSourceMatch(
                new File("src/dr/evomodel/treedatalikelihood/continuous"),
                WILDCARD_IMPORT,
                "OUCanonical");
        assertNoSourceMatch(
                new File("src/dr/evomodel/continuous/canonical"),
                WILDCARD_IMPORT);
        assertNoSourceMatch(
                new File("src/dr/evomodel/continuous/ou/canonical"),
                WILDCARD_IMPORT);
        assertNoSourceMatch(
                new File("src/dr/evomodel/continuous/ou/blockdiagonal"),
                WILDCARD_IMPORT);
        assertNoSourceMatch(
                new File("src/dr/evomodel/treedatalikelihood/hmc"),
                WILDCARD_IMPORT,
                "Canonical");
    }

    public void testCanonicalSubpackagesRespectArchitectureBoundaries() throws IOException {
        assertNoImportContaining(
                new File("src/dr/evomodel/treedatalikelihood/continuousNew/canonical/traversal"),
                "dr.evomodel.treedatalikelihood.continuousNew.canonical.gradient");
        assertNoImportContaining(
                new File("src/dr/evomodel/treedatalikelihood/continuousNew/canonical/traversal"),
                "dr.evomodel.treedatalikelihood.continuousNew.canonical.adapter");
        assertNoImportContaining(
                new File("src/dr/evomodel/treedatalikelihood/continuousNew/canonical/traversal"),
                BLOCK_DIAGONAL_BACKEND_PACKAGE);

        assertNoImportContaining(
                new File("src/dr/evomodel/treedatalikelihood/continuousNew/canonical/message"),
                "dr.evolution.tree");
        assertNoImportContaining(
                new File("src/dr/evomodel/treedatalikelihood/continuousNew/canonical/message"),
                "dr.evomodel.treedatalikelihood.continuousNew.canonical.traversal");
        assertNoImportContaining(
                new File("src/dr/evomodel/treedatalikelihood/continuousNew/canonical/message"),
                "dr.evomodel.treedatalikelihood.continuousNew.canonical.adapter");

        assertNoImportContaining(
                new File("src/dr/evomodel/treedatalikelihood/continuousNew/canonical/gradient"),
                "dr.evomodel.treedatalikelihood.continuousNew.canonical.adapter");
        assertNoImportContaining(
                new File("src/dr/evomodel/treedatalikelihood/continuousNew/canonical/gradient"),
                BLOCK_DIAGONAL_BACKEND_PACKAGE);

        assertNoImportContaining(
                new File("src/dr/evomodel/continuous/ou/canonical"),
                "dr.evomodel.treedatalikelihood.continuousNew.canonical.adapter");
        assertNoImportContaining(
                new File("src/dr/evomodel/continuous/ou/canonical"),
                "dr.evomodel.treedatalikelihood.continuousNew.canonical.traversal");
        assertNoImportContaining(
                new File("src/dr/evomodel/continuous/ou/canonical"),
                BLOCK_DIAGONAL_BACKEND_PACKAGE);

        assertNoImportContaining(
                new File("src/dr/evomodel/continuous/ou/blockdiagonal"),
                "dr.evomodel.treedatalikelihood.continuousNew.canonical.adapter");
        assertNoImportContaining(
                new File("src/dr/evomodel/continuous/ou/blockdiagonal"),
                "dr.evomodel.treedatalikelihood.continuousNew.canonical.traversal");
        assertNoImportContaining(
                new File("src/dr/evomodel/continuous/ou/blockdiagonal"),
                "dr.evomodel.treedatalikelihood.continuousNew.canonical.gradient");

        // dr.evomodel.continuous.canonical (math, contribution, workspace, scheduling,
        // BMProcessModel/ContinuousProcessModel) is the relocated pure-math half of the former
        // continuousNew.canonical.* subpackages (see private-code/docs/
        // continuousNew_tree_math_separation_plan.md, "Dependency Rules"): it must never depend
        // back on the tree-dependent traversal/adapter code it was split out of, or on the tree
        // itself.
        assertNoImportContaining(
                new File("src/dr/evomodel/continuous/canonical"),
                "dr.evomodel.treedatalikelihood.continuousNew.canonical.traversal");
        assertNoImportContaining(
                new File("src/dr/evomodel/continuous/canonical"),
                "dr.evomodel.treedatalikelihood.continuousNew.canonical.adapter");
        assertNoImportContaining(
                new File("src/dr/evomodel/continuous/canonical"),
                "dr.evolution.tree");
    }

    public void testCanonicalTreeHotPathsUseFlatMatrixOps() throws IOException {
        assertFalse(
                "Legacy dense matrix compatibility wrapper should not exist",
                new File(
                        "src/dr/evomodel/treedatalikelihood/continuousNew/canonical/message",
                        "Gaussian" + "MatrixOps.java")
                        .exists());
    }

    public void testCanonicalGradientBackendDoesNotDependOnLegacyBranchGradientPath() throws IOException {
        assertNoImportContaining(
                new File("src/dr/evomodel/treedatalikelihood/continuousNew/canonical/adapter"),
                "dr.evomodel.treedatalikelihood.continuous.BranchSpecificGradient");
        assertNoImportContaining(
                new File("src/dr/evomodel/treedatalikelihood/continuousNew/canonical/adapter"),
                "dr.evomodel.treedatalikelihood.continuous.ContinuousTraitGradientForBranch");
        assertNoImportContaining(
                new File("src/dr/evomodel/treedatalikelihood/continuousNew/canonical/gradient"),
                "dr.evomodel.treedatalikelihood.continuous.BranchSpecificGradient");
        assertNoImportContaining(
                new File("src/dr/evomodel/treedatalikelihood/continuousNew/canonical/gradient"),
                "dr.evomodel.treedatalikelihood.continuous.ContinuousTraitGradientForBranch");
    }

    public void testCanonicalGradientInternalsAvoidDoubleArrayMatrices() throws IOException {
        assertNoSourceMatch(
                new File("src/dr/evomodel/treedatalikelihood/continuousNew/canonical/gradient"),
                "double[][]");
        assertNoSourceMatch(
                new File("src/dr/evomodel/continuous/ou/canonical"),
                "double[][]");
    }

    public void testBlockDiagonalHotPathUsesFillPullbacks() throws IOException {
        assertNoSourceMatch(
                new File("src/dr/evomodel/continuous/ou/blockdiagonal"),
                "pullBackGradientFlat(");
        assertNoSourceMatch(
                new File("src/dr/evomodel/treedatalikelihood/continuousNew/canonical/gradient/SpecializedCanonicalSelectionGradientPullback.java"),
                "projectDenseGradient(");
    }

    public void testBlockDiagonalBackendAllocatesDoubleBuffersOnlyDuringConstruction() throws IOException {
        assertNoMethodLocalDoubleArrayAllocation(
                new File("src/dr/evomodel/continuous/ou/blockdiagonal"));
    }

    private static void assertNoSourceMatch(final File root,
                                            final String forbidden) throws IOException {
        assertNoSourceMatch(root, forbidden, null);
    }

    private static void assertNoSourceMatch(final File root,
                                            final String forbidden,
                                            final String fileNamePrefix) throws IOException {
        if (!root.exists()) {
            fail("Source root does not exist: " + root.getPath());
        }
        scan(root, forbidden, fileNamePrefix);
    }

    private static void scan(final File file,
                             final String forbidden,
                             final String fileNamePrefix) throws IOException {
        if (file.isDirectory()) {
            final File[] children = file.listFiles();
            if (children == null) {
                return;
            }
            for (File child : children) {
                scan(child, forbidden, fileNamePrefix);
            }
            return;
        }
        if (!file.getName().endsWith(".java")) {
            return;
        }
        if (fileNamePrefix != null && !file.getName().startsWith(fileNamePrefix)) {
            return;
        }
        final BufferedReader reader = new BufferedReader(new FileReader(file));
        try {
            String line;
            int lineNumber = 0;
            while ((line = reader.readLine()) != null) {
                lineNumber++;
                if (matches(line, forbidden)) {
                    fail(file.getPath() + ":" + lineNumber + " contains forbidden text: " + forbidden);
                }
            }
        } finally {
            reader.close();
        }
    }

    private static void assertNoImportContaining(final File root,
                                                 final String forbiddenImport) throws IOException {
        if (!root.exists()) {
            fail("Source root does not exist: " + root.getPath());
        }
        scanImports(root, forbiddenImport);
    }

    private static void scanImports(final File file,
                                    final String forbiddenImport) throws IOException {
        if (file.isDirectory()) {
            final File[] children = file.listFiles();
            if (children == null) {
                return;
            }
            for (File child : children) {
                scanImports(child, forbiddenImport);
            }
            return;
        }
        if (!file.getName().endsWith(".java")) {
            return;
        }
        final BufferedReader reader = new BufferedReader(new FileReader(file));
        try {
            String line;
            int lineNumber = 0;
            while ((line = reader.readLine()) != null) {
                lineNumber++;
                final String trimmed = line.trim();
                if (trimmed.startsWith("import ") && trimmed.contains(forbiddenImport)) {
                    fail(file.getPath() + ":" + lineNumber + " imports forbidden boundary: " + forbiddenImport);
                }
            }
        } finally {
            reader.close();
        }
    }

    private static boolean matches(final String line, final String forbidden) {
        if (WILDCARD_IMPORT.equals(forbidden)) {
            final String trimmed = line.trim();
            return trimmed.startsWith("import ") && trimmed.endsWith(".*;");
        }
        if (SYSTEM_PRINT.equals(forbidden)) {
            return line.contains("System.err") || line.contains("System.out");
        }
        return line.contains(forbidden);
    }

    private static void assertNoMethodLocalDoubleArrayAllocation(final File root) throws IOException {
        if (!root.exists()) {
            fail("Source root does not exist: " + root.getPath());
        }
        scanDoubleArrayAllocations(root);
    }

    private static void scanDoubleArrayAllocations(final File file) throws IOException {
        if (file.isDirectory()) {
            final File[] children = file.listFiles();
            if (children == null) {
                return;
            }
            for (File child : children) {
                scanDoubleArrayAllocations(child);
            }
            return;
        }
        if (!file.getName().endsWith(".java")) {
            return;
        }
        final BufferedReader reader = new BufferedReader(new FileReader(file));
        try {
            String line;
            int lineNumber = 0;
            // Java forbids bare executable statements directly inside a class/interface/enum
            // body - anything declared one level inside such a body is necessarily a field,
            // method, or nested type, never a hot-path local variable. So track, brace by brace,
            // whether we are currently directly inside a type body ("true") or inside a
            // method/constructor/block body ("false"); allocations are only flagged when the
            // innermost context is a method body.
            final java.util.Deque<Boolean> contextStack = new java.util.ArrayDeque<>();
            // One-time static array-builder methods (e.g. a "private static double[]
            // buildInverseFactorials(...)" called exactly once to seed a static final constant
            // table) are the method-body equivalent of a constructor-owned allocation; recognized
            // by their signature line and suppressed until the matching closing brace.
            boolean insideStaticArrayFactory = false;
            int staticArrayFactoryDepth = -1;
            while ((line = reader.readLine()) != null) {
                lineNumber++;
                final String trimmed = line.trim();

                if (trimmed.startsWith("//") || trimmed.isEmpty()) {
                    continue;
                }

                final boolean atTypeBodyContext = contextStack.isEmpty() || contextStack.peek();
                // A "this.field = new double[...]" assignment is a constructor-owned allocation
                // even though the constructor body itself is a method-like (non-type) context.
                final boolean isConstructorFieldAssignment = trimmed.startsWith("this.");

                if (!insideStaticArrayFactory
                        && !atTypeBodyContext
                        && !isConstructorFieldAssignment
                        && trimmed.contains("new double[")) {
                    fail(file.getPath() + ":" + lineNumber
                            + " allocates a double[] outside constructor-owned workspace fields");
                }

                final boolean opensType = trimmed.matches(
                        "(private |protected |public |static |final |abstract )*"
                                + "(class|interface|enum) \\w+.*\\{\\s*");
                final boolean opensBrace = trimmed.endsWith("{");
                final boolean closesBrace = trimmed.equals("}") || trimmed.equals("};");

                if (insideStaticArrayFactory && closesBrace && contextStack.size() == staticArrayFactoryDepth) {
                    insideStaticArrayFactory = false;
                    staticArrayFactoryDepth = -1;
                }

                if (opensBrace) {
                    contextStack.push(opensType);
                    if (!opensType && trimmed.matches("private static \\S*\\[\\]\\s+\\w+\\(.*\\)\\s*\\{")) {
                        insideStaticArrayFactory = true;
                        staticArrayFactoryDepth = contextStack.size();
                    }
                } else if (closesBrace && !contextStack.isEmpty()) {
                    contextStack.pop();
                }
            }
        } finally {
            reader.close();
        }
    }

    private static String readFile(final File file) throws IOException {
        final StringBuilder builder = new StringBuilder();
        final BufferedReader reader = new BufferedReader(new FileReader(file));
        try {
            String line;
            while ((line = reader.readLine()) != null) {
                builder.append(line).append('\n');
            }
        } finally {
            reader.close();
        }
        return builder.toString();
    }
}

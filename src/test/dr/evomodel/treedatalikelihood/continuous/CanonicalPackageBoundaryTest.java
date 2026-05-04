package test.dr.evomodel.treedatalikelihood.continuous;

import junit.framework.TestCase;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

public final class CanonicalPackageBoundaryTest extends TestCase {

    private static final String ORTHOGONAL_BACKEND_PACKAGE =
            "dr.evomodel.continuous.ou.orthogonalblockdiagonal";
    private static final String LEGACY_GAUSSIAN_MATRIX_OPS =
            "dr.evomodel.treedatalikelihood.continuous.canonical.message.GaussianMatrixOps";
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

    public void testCanonicalArchitectureNoteLocksOrthogonalBlockPullbackContract() throws IOException {
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

    public void testCanonicalTreeLayerDoesNotImportOrthogonalBackend() throws IOException {
        assertNoSourceMatch(
                new File("src/dr/evomodel/treedatalikelihood/continuous/canonical"),
                ORTHOGONAL_BACKEND_PACKAGE);
        assertNoSourceMatch(
                new File("src/dr/evomodel/treedatalikelihood/continuous"),
                ORTHOGONAL_BACKEND_PACKAGE,
                "OUCanonical");
    }

    public void testCanonicalMachineryDoesNotUseDirectSystemPrints() throws IOException {
        assertNoSourceMatch(
                new File("src/dr/evomodel/treedatalikelihood/continuous/canonical"),
                SYSTEM_PRINT);
        assertNoSourceMatch(
                new File("src/dr/evomodel/continuous/ou/canonical"),
                SYSTEM_PRINT);
        assertNoSourceMatch(
                new File("src/dr/evomodel/continuous/ou/orthogonalblockdiagonal"),
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
                new File("src/dr/evomodel/treedatalikelihood/continuous/canonical"),
                WILDCARD_IMPORT);
        assertNoSourceMatch(
                new File("src/dr/evomodel/treedatalikelihood/continuous"),
                WILDCARD_IMPORT,
                "OUCanonical");
        assertNoSourceMatch(
                new File("src/dr/evomodel/continuous/ou/canonical"),
                WILDCARD_IMPORT);
        assertNoSourceMatch(
                new File("src/dr/evomodel/continuous/ou/orthogonalblockdiagonal"),
                WILDCARD_IMPORT);
        assertNoSourceMatch(
                new File("src/dr/evomodel/treedatalikelihood/hmc"),
                WILDCARD_IMPORT,
                "Canonical");
    }

    public void testCanonicalSubpackagesRespectArchitectureBoundaries() throws IOException {
        assertNoImportContaining(
                new File("src/dr/evomodel/treedatalikelihood/continuous/canonical/traversal"),
                "dr.evomodel.treedatalikelihood.continuous.canonical.gradient");
        assertNoImportContaining(
                new File("src/dr/evomodel/treedatalikelihood/continuous/canonical/traversal"),
                "dr.evomodel.treedatalikelihood.continuous.canonical.adapter");
        assertNoImportContaining(
                new File("src/dr/evomodel/treedatalikelihood/continuous/canonical/traversal"),
                "dr.evomodel.continuous.ou.orthogonalblockdiagonal");

        assertNoImportContaining(
                new File("src/dr/evomodel/treedatalikelihood/continuous/canonical/message"),
                "dr.evolution.tree");
        assertNoImportContaining(
                new File("src/dr/evomodel/treedatalikelihood/continuous/canonical/message"),
                "dr.evomodel.treedatalikelihood.continuous.canonical.traversal");
        assertNoImportContaining(
                new File("src/dr/evomodel/treedatalikelihood/continuous/canonical/message"),
                "dr.evomodel.treedatalikelihood.continuous.canonical.adapter");

        assertNoImportContaining(
                new File("src/dr/evomodel/treedatalikelihood/continuous/canonical/gradient"),
                "dr.evomodel.treedatalikelihood.continuous.canonical.adapter");
        assertNoImportContaining(
                new File("src/dr/evomodel/treedatalikelihood/continuous/canonical/gradient"),
                "dr.evomodel.continuous.ou.orthogonalblockdiagonal");

        assertNoImportContaining(
                new File("src/dr/evomodel/continuous/ou/canonical"),
                "dr.evomodel.treedatalikelihood.continuous.canonical.adapter");
        assertNoImportContaining(
                new File("src/dr/evomodel/continuous/ou/canonical"),
                "dr.evomodel.treedatalikelihood.continuous.canonical.traversal");
        assertNoImportContaining(
                new File("src/dr/evomodel/continuous/ou/canonical"),
                "dr.evomodel.continuous.ou.orthogonalblockdiagonal");

        assertNoImportContaining(
                new File("src/dr/evomodel/continuous/ou/orthogonalblockdiagonal"),
                "dr.evomodel.treedatalikelihood.continuous.canonical.adapter");
        assertNoImportContaining(
                new File("src/dr/evomodel/continuous/ou/orthogonalblockdiagonal"),
                "dr.evomodel.treedatalikelihood.continuous.canonical.traversal");
        assertNoImportContaining(
                new File("src/dr/evomodel/continuous/ou/orthogonalblockdiagonal"),
                "dr.evomodel.treedatalikelihood.continuous.canonical.gradient");
    }

    public void testCanonicalTreeHotPathsUseFlatMatrixOps() throws IOException {
        assertNoImportContaining(
                new File("src/dr/evomodel/treedatalikelihood/continuous/canonical/adapter"),
                LEGACY_GAUSSIAN_MATRIX_OPS);
        assertNoImportContaining(
                new File("src/dr/evomodel/treedatalikelihood/continuous/canonical/contribution"),
                LEGACY_GAUSSIAN_MATRIX_OPS);
        assertNoImportContaining(
                new File("src/dr/evomodel/treedatalikelihood/continuous/canonical/gradient"),
                LEGACY_GAUSSIAN_MATRIX_OPS);
        assertNoImportContaining(
                new File("src/dr/evomodel/treedatalikelihood/continuous/canonical/traversal"),
                LEGACY_GAUSSIAN_MATRIX_OPS);
        assertNoImportContaining(
                new File("src/dr/evomodel/treedatalikelihood/continuous/canonical/workspace"),
                LEGACY_GAUSSIAN_MATRIX_OPS);
        assertNoImportContaining(
                new File("src/dr/evomodel/continuous/ou/canonical"),
                LEGACY_GAUSSIAN_MATRIX_OPS);
        assertNoImportContaining(
                new File("src/dr/evomodel/continuous/ou/orthogonalblockdiagonal"),
                LEGACY_GAUSSIAN_MATRIX_OPS);
    }

    public void testCanonicalGradientBackendDoesNotDependOnLegacyBranchGradientPath() throws IOException {
        assertNoImportContaining(
                new File("src/dr/evomodel/treedatalikelihood/continuous/canonical/adapter"),
                "dr.evomodel.treedatalikelihood.continuous.BranchSpecificGradient");
        assertNoImportContaining(
                new File("src/dr/evomodel/treedatalikelihood/continuous/canonical/adapter"),
                "dr.evomodel.treedatalikelihood.continuous.ContinuousTraitGradientForBranch");
        assertNoImportContaining(
                new File("src/dr/evomodel/treedatalikelihood/continuous/canonical/gradient"),
                "dr.evomodel.treedatalikelihood.continuous.BranchSpecificGradient");
        assertNoImportContaining(
                new File("src/dr/evomodel/treedatalikelihood/continuous/canonical/gradient"),
                "dr.evomodel.treedatalikelihood.continuous.ContinuousTraitGradientForBranch");
    }

    public void testCanonicalGradientInternalsAvoidDoubleArrayMatrices() throws IOException {
        assertNoSourceMatch(
                new File("src/dr/evomodel/treedatalikelihood/continuous/canonical/gradient"),
                "double[][]");
        assertNoSourceMatch(
                new File("src/dr/evomodel/continuous/ou/canonical"),
                "double[][]");
    }

    public void testOrthogonalBlockHotPathUsesFillPullbacks() throws IOException {
        assertNoSourceMatch(
                new File("src/dr/evomodel/continuous/ou/orthogonalblockdiagonal"),
                "pullBackGradientFlat(");
        assertNoSourceMatch(
                new File("src/dr/evomodel/treedatalikelihood/continuous/canonical/gradient/SpecializedCanonicalSelectionGradientPullback.java"),
                "projectDenseGradient(");
    }

    public void testOrthogonalBlockBackendAllocatesDoubleBuffersOnlyDuringConstruction() throws IOException {
        assertNoMethodLocalDoubleArrayAllocation(
                new File("src/dr/evomodel/continuous/ou/orthogonalblockdiagonal"));
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
            while ((line = reader.readLine()) != null) {
                lineNumber++;
                final String trimmed = line.trim();
                if (trimmed.contains("new double[") && !trimmed.startsWith("this.")) {
                    fail(file.getPath() + ":" + lineNumber
                            + " allocates a double[] outside constructor-owned workspace fields");
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

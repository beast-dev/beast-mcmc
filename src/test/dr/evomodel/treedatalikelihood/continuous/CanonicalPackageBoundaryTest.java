package test.dr.evomodel.treedatalikelihood.continuous;

import junit.framework.TestCase;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

public final class CanonicalPackageBoundaryTest extends TestCase {

    private static final String ORTHOGONAL_BACKEND_PACKAGE =
            "dr.evomodel.continuous.ou.orthogonalblockdiagonal";
    private static final String WILDCARD_IMPORT = "wildcard import";

    public void testCanonicalTreeLayerDoesNotImportOrthogonalBackend() throws IOException {
        assertNoSourceMatch(
                new File("src/dr/evomodel/treedatalikelihood/continuous/canonical"),
                ORTHOGONAL_BACKEND_PACKAGE);
    }

    public void testCanonicalMachineryDoesNotUseWildcardImports() throws IOException {
        assertNoSourceMatch(
                new File("src/dr/evomodel/treedatalikelihood/continuous/canonical"),
                WILDCARD_IMPORT);
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

    private static boolean matches(final String line, final String forbidden) {
        if (WILDCARD_IMPORT.equals(forbidden)) {
            final String trimmed = line.trim();
            return trimmed.startsWith("import ") && trimmed.endsWith(".*;");
        }
        return line.contains(forbidden);
    }
}

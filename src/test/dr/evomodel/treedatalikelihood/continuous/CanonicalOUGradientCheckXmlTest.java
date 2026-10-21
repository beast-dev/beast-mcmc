package test.dr.evomodel.treedatalikelihood.continuous;

import dr.app.beast.BeastParser;
import dr.evomodel.treedatalikelihood.continuousNew.CanonicalOUBranchRateGradient;
import dr.evomodel.treedatalikelihood.continuousNew.CanonicalOUNodeHeightGradient;
import dr.inference.mcmc.MCMC;
import dr.inference.operators.hmc.HamiltonianMonteCarloOperator;
import dr.inference.operators.hmc.ReflectiveHamiltonianMonteCarloOperator;
import dr.xml.XMLObject;
import junit.framework.TestCase;

import java.io.File;
import java.io.FileReader;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.Map;

public final class CanonicalOUGradientCheckXmlTest extends TestCase {

    private static final String BRANCH_RATE_XML =
            "tests/TestXML/testCanonicalOUBranchRateGradientCheck.xml";
    private static final String NODE_HEIGHT_XML =
            "tests/TestXML/testCanonicalOUNodeHeightGradientCheck.xml";
    private static final String BLOCK_BASIS_PARTIAL_XML =
            "tests/TestXML/testCanonicalOUBlockBasisPartialTips.xml";
    private static final String BLOCK_BASIS_ANOLIS_FULL_TIPS_XML =
            "tests/TestXML/testCanonicalOUBlockBasisAnolisFullTips.xml";

    public void testStrictClockBranchRateXmlRunsHmcGradientCheck() throws Exception {
        final ParsedXml parsed = parseXml(BRANCH_RATE_XML);

        assertNativeObject(parsed.objects, "canonical.clock.gradient", CanonicalOUBranchRateGradient.class);
        assertNativeObject(parsed.objects, "canonicalClockHMC", HamiltonianMonteCarloOperator.class);
        assertGradientCheckAttributes(parsed.file, "100", "1e-4");

        runMcmc(parsed.objects);
    }

    public void testNodeHeightXmlRunsReflectiveHmcGradientCheck() throws Exception {
        final ParsedXml parsed = parseXml(NODE_HEIGHT_XML);

        assertNativeObject(parsed.objects, "canonical.nodeHeight.gradient", CanonicalOUNodeHeightGradient.class);
        assertNativeObject(parsed.objects, "canonicalNodeHeightHMC", ReflectiveHamiltonianMonteCarloOperator.class);
        assertGradientCheckAttributes(parsed.file, "100", "1e-4");

        runMcmc(parsed.objects);
    }

    public void testBlockBasisPartialTipXmlRunsBranchRateGradientCheck() throws Exception {
        final ParsedXml parsed = parseXml(BLOCK_BASIS_PARTIAL_XML);

        assertNativeObject(parsed.objects, "block.clock.gradient", CanonicalOUBranchRateGradient.class);
        assertNativeObject(parsed.objects, "blockClockHMC", HamiltonianMonteCarloOperator.class);
        assertGradientCheckAttributes(parsed.file, "100", "1e-4");

        runMcmc(parsed.objects);
    }

    public void testBlockBasisAnolisFullTipXmlRunsHmcGradientCheck() throws Exception {
        final ParsedXml parsed = parseXml(BLOCK_BASIS_ANOLIS_FULL_TIPS_XML);

        assertNativeObject(parsed.objects, "anoleBlockBasisFullTipHMC", HamiltonianMonteCarloOperator.class);
        assertGradientCheckAttributes(parsed.file, "10", "1e-3");

        runMcmc(parsed.objects);
    }

    private static ParsedXml parseXml(final String path) throws Exception {
        final File file = new File(path);
        assertTrue("missing TestXML file: " + file.getAbsolutePath(), file.isFile());

        final BeastParser parser = new BeastParser(new String[0], null, false, true, true, null);
        final Map<String, XMLObject> objects = parser.parse(new FileReader(file), true);
        return new ParsedXml(file, objects);
    }

    private static void assertNativeObject(final Map<String, XMLObject> objects,
                                           final String id,
                                           final Class<?> expectedType) {
        final XMLObject object = objects.get(id);
        assertNotNull("XML object should be stored by id: " + id, object);
        assertTrue("XML object should be " + expectedType.getSimpleName() + ": " + id,
                expectedType.isInstance(object.getNativeObject()));
    }

    private static void assertGradientCheckAttributes(final File file,
                                                      final String expectedCount,
                                                      final String expectedTolerance) throws Exception {
        final String xml = new String(Files.readAllBytes(file.toPath()), StandardCharsets.UTF_8);
        assertTrue("HMC XML should use the stress gradient check count",
                xml.contains("gradientCheckCount=\"" + expectedCount + "\""));
        assertTrue("HMC XML should use the stress gradient check tolerance",
                xml.contains("gradientCheckTolerance=\"" + expectedTolerance + "\""));
    }

    private static void runMcmc(final Map<String, XMLObject> objects) {
        final XMLObject object = objects.get("mcmc");
        assertNotNull("XML object should be stored by id: mcmc", object);
        assertTrue("mcmc object should parse to MCMC", object.getNativeObject() instanceof MCMC);
        ((MCMC) object.getNativeObject()).chain();
    }

    private static final class ParsedXml {
        final File file;
        final Map<String, XMLObject> objects;

        private ParsedXml(final File file,
                          final Map<String, XMLObject> objects) {
            this.file = file;
            this.objects = objects;
        }
    }
}

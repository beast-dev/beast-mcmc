package test.dr.evomodel.treedatalikelihood.continuous;

import dr.app.beast.BeastParser;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.CanonicalOUBranchRateGradient;
import dr.evomodel.treedatalikelihood.continuous.CanonicalOUNodeHeightGradient;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.xml.XMLObject;
import junit.framework.TestCase;

import java.io.StringReader;
import java.lang.reflect.Field;
import java.util.Map;

public final class ContinuousDataLikelihoodParserSelectionTest extends TestCase {

    private static final String CDI_BM_ENGINE =
            "dr.evomodel.treedatalikelihood.continuous.CdiBMContinuousTreeLikelihoodEngine";
    private static final String MOMENT_BM_ENGINE =
            "dr.evomodel.treedatalikelihood.continuous.MomentBMContinuousTreeLikelihoodEngine";
    private static final String CANONICAL_OU_ENGINE =
            "dr.evomodel.treedatalikelihood.continuous.CanonicalOUContinuousTreeLikelihoodEngine";

    public void testParserCanonicalCompatibilityBMSelectsCdiMomentBackend() throws Exception {
        final ContinuousDataLikelihoodDelegate delegate = parseDelegate(bmXml("canonical"));

        assertTrue("legacy canonical BM spelling should select CDI BM", delegate.usesCdiBMLikelihood());
        assertFalse("BM should not select canonical OU", delegate.usesCanonicalOULikelihood());
        assertSelectedEngine(delegate, CDI_BM_ENGINE);
    }

    public void testParserExplicitCdiMomentBMSelectsCdiMomentBackend() throws Exception {
        final ContinuousDataLikelihoodDelegate delegate = parseDelegate(bmXml("cdi-moment"));

        assertTrue("explicit CDI/moment BM spelling should select CDI BM", delegate.usesCdiBMLikelihood());
        assertFalse("explicit CDI/moment BM should not select native BM", delegate.usesMomentBMLikelihood());
        assertFalse("BM should not select canonical OU", delegate.usesCanonicalOULikelihood());
        assertSelectedEngine(delegate, CDI_BM_ENGINE);
    }

    public void testParserExplicitNativeMomentBMSelectsMomentBackendAndMatchesCdi() throws Exception {
        final TreeDataLikelihood cdiLikelihood = parseLikelihood(bmXml("cdi-moment"));
        final TreeDataLikelihood momentLikelihood = parseLikelihood(bmXml("bm-moment"));
        final ContinuousDataLikelihoodDelegate momentDelegate =
                (ContinuousDataLikelihoodDelegate) momentLikelihood.getDataLikelihoodDelegate();

        assertTrue("explicit native moment BM should select native BM", momentDelegate.usesMomentBMLikelihood());
        assertFalse("explicit native moment BM should not select CDI BM", momentDelegate.usesCdiBMLikelihood());
        assertFalse("BM should not select canonical OU", momentDelegate.usesCanonicalOULikelihood());
        assertSelectedEngine(momentDelegate, MOMENT_BM_ENGINE);
        assertEquals("native moment BM should match CDI/moment BM",
                cdiLikelihood.getLogLikelihood(), momentLikelihood.getLogLikelihood(), 1.0e-10);
    }

    public void testParserCanonicalCompatibilityOUSelectsCanonicalMessageBackend() throws Exception {
        final ContinuousDataLikelihoodDelegate delegate = parseDelegate(ouXml("canonical"));

        assertTrue("legacy canonical OU spelling should select canonical OU", delegate.usesCanonicalOULikelihood());
        assertFalse("OU should not select CDI BM", delegate.usesCdiBMLikelihood());
        assertSelectedEngine(delegate, CANONICAL_OU_ENGINE);
    }

    public void testParserExplicitCanonicalMessageOUSelectsCanonicalMessageBackend() throws Exception {
        final ContinuousDataLikelihoodDelegate delegate = parseDelegate(ouXml("canonical-message"));

        assertTrue("explicit canonical-message OU spelling should select canonical OU", delegate.usesCanonicalOULikelihood());
        assertFalse("OU should not select CDI BM", delegate.usesCdiBMLikelihood());
        assertSelectedEngine(delegate, CANONICAL_OU_ENGINE);
    }

    public void testParserCanonicalMessageOUBranchProcessesSelectsCanonicalMessageBackend() throws Exception {
        final ContinuousDataLikelihoodDelegate delegate = parseDelegate(heterogeneousOuXml("canonical-message"));

        assertTrue("heterogeneous OU branch processes should select canonical OU",
                delegate.usesCanonicalOULikelihood());
        assertFalse("heterogeneous OU branch processes should not select CDI BM",
                delegate.usesCdiBMLikelihood());
        assertSelectedEngine(delegate, CANONICAL_OU_ENGINE);
    }

    public void testParserRejectsCanonicalMessageForBM() {
        assertParseFails(
                bmXml("canonical-message"),
                "canonical-message OU backend and supports only OU models");
    }

    public void testParserRejectsCdiMomentForOU() {
        assertParseFails(
                ouXml("cdi-moment"),
                "CDI/moment-form BM production backend and supports only homogeneous BM");
    }

    public void testParserRejectsNativeMomentBMForOU() {
        assertParseFails(
                ouXml("bm-moment"),
                "native moment-form BM backend and supports only homogeneous BM");
    }

    public void testParserRejectsNaturalCanonicalBMProductionSelection() {
        assertParseFails(
                bmXml("natural-canonical-bm"),
                "experimental natural-parameter BM backend");
    }

    public void testParserRejectsOUBranchProcessesForLegacyBackend() {
        assertParseFails(
                heterogeneousOuXml("legacy"),
                "requires implementation=\"canonical-message\" or implementation=\"canonical\"");
    }

    public void testParserCanonicalOUBranchRateGradientSelectsCanonicalProvider() throws Exception {
        final Object gradient = parseNativeObject(
                withGradientElement(
                        ouXml("canonical-message"),
                        "    <branchRateGradient id=\"clockGradient\" traitName=\"X\">"
                                + "        <traitDataLikelihood idref=\"likelihood\"/>"
                                + "    </branchRateGradient>"),
                "clockGradient");

        assertTrue("canonical OU branchRateGradient should select canonical provider",
                gradient instanceof CanonicalOUBranchRateGradient);
    }

    public void testParserCanonicalOUNodeHeightGradientSelectsCanonicalProvider() throws Exception {
        final Object gradient = parseNativeObject(
                withGradientElement(
                        ouXml("canonical-message"),
                        "    <nodeHeightGradient id=\"heightGradient\" traitName=\"X\">"
                                + "        <traitDataLikelihood idref=\"likelihood\"/>"
                                + "    </nodeHeightGradient>"),
                "heightGradient");

        assertTrue("canonical OU nodeHeightGradient should select canonical provider",
                gradient instanceof CanonicalOUNodeHeightGradient);
    }

    public void testParserRejectsCanonicalOUBranchRateGradientWithoutRateParameter() {
        assertNativeObjectParseFails(
                withGradientElement(
                        ouXmlWithoutTraitRate("canonical-message"),
                        "    <branchRateGradient id=\"clockGradient\" traitName=\"X\">"
                                + "        <traitDataLikelihood idref=\"likelihood\"/>"
                                + "    </branchRateGradient>"),
                "clockGradient",
                "explicit rate parameter");
    }

    public void testParserRejectsCanonicalOUNodeHeightGradientForTreeDependentRates() {
        assertNativeObjectParseFails(
                withGradientElement(
                        ouXmlWithTreeTimeScaledTraitRate("canonical-message"),
                        "    <nodeHeightGradient id=\"heightGradient\" traitName=\"X\">"
                                + "        <traitDataLikelihood idref=\"likelihood\"/>"
                                + "    </nodeHeightGradient>"),
                "heightGradient",
                "branch rates that do not depend on node heights");
    }

    private static ContinuousDataLikelihoodDelegate parseDelegate(final String xml) throws Exception {
        final TreeDataLikelihood likelihood = parseLikelihood(xml);
        return (ContinuousDataLikelihoodDelegate) likelihood.getDataLikelihoodDelegate();
    }

    private static TreeDataLikelihood parseLikelihood(final String xml) throws Exception {
        return (TreeDataLikelihood) parseNativeObject(xml, "likelihood");
    }

    private static Object parseNativeObject(final String xml, final String id) throws Exception {
        final BeastParser parser = new BeastParser(new String[0], null, false, true, true, null);
        final Map<String, XMLObject> objects = parser.parse(new StringReader(xml), true);
        final XMLObject object = objects.get(id);
        assertNotNull(id + " object should be stored by id", object);
        return object.getNativeObject();
    }

    private static void assertParseFails(final String xml, final String expectedMessage) {
        try {
            parseDelegate(xml);
            fail("Expected parser rejection containing: " + expectedMessage);
        } catch (Exception expected) {
            assertTrue("Expected parser message to contain '" + expectedMessage + "' but got: "
                            + messageChain(expected),
                    messageChain(expected).contains(expectedMessage));
        }
    }

    private static void assertNativeObjectParseFails(final String xml,
                                                     final String id,
                                                     final String expectedMessage) {
        try {
            parseNativeObject(xml, id);
            fail("Expected parser rejection containing: " + expectedMessage);
        } catch (Exception expected) {
            assertTrue("Expected parser message to contain '" + expectedMessage + "' but got: "
                            + messageChain(expected),
                    messageChain(expected).contains(expectedMessage));
        }
    }

    private static void assertSelectedEngine(final ContinuousDataLikelihoodDelegate delegate,
                                             final String expectedClassName) throws Exception {
        final Object engine = privateFieldValue(delegate, "likelihoodEngine");
        assertNotNull("likelihood engine should be selected", engine);
        assertEquals(expectedClassName, engine.getClass().getName());
    }

    private static Object privateFieldValue(final Object target, final String fieldName) throws Exception {
        final Field field = declaredField(target.getClass(), fieldName);
        field.setAccessible(true);
        return field.get(target);
    }

    private static Field declaredField(final Class<?> initialType, final String fieldName) throws NoSuchFieldException {
        Class<?> type = initialType;
        while (type != null) {
            try {
                return type.getDeclaredField(fieldName);
            } catch (NoSuchFieldException ignored) {
                type = type.getSuperclass();
            }
        }
        throw new NoSuchFieldException(fieldName);
    }

    private static String messageChain(final Throwable throwable) {
        final StringBuilder builder = new StringBuilder();
        Throwable current = throwable;
        while (current != null) {
            if (current.getMessage() != null) {
                if (builder.length() > 0) {
                    builder.append(" | ");
                }
                builder.append(current.getMessage());
            }
            current = current.getCause();
        }
        return builder.toString();
    }

    private static String bmXml(final String implementation) {
        return xml(implementation, "");
    }

    private static String ouXml(final String implementation) {
        return xml(implementation,
                "        <optimalTraits>"
                        + "            <strictClockBranchRates><rate><parameter id=\"opt.1\" value=\"1.0\"/></rate></strictClockBranchRates>"
                        + "            <strictClockBranchRates><rate><parameter id=\"opt.2\" value=\"2.0\"/></rate></strictClockBranchRates>"
                        + "        </optimalTraits>"
                        + "        <strengthOfSelectionMatrix>"
                        + "            <diagonalMatrix id=\"selection\">"
                        + "                <parameter id=\"alpha.diag\" value=\"1.0 1.5\" lower=\"0.0 0.0\"/>"
                        + "            </diagonalMatrix>"
                        + "        </strengthOfSelectionMatrix>");
    }

    private static String ouXmlWithoutTraitRate(final String implementation) {
        return xml(implementation,
                "        <optimalTraits>"
                        + "            <strictClockBranchRates><rate><parameter id=\"no.rate.opt.1\" value=\"1.0\"/></rate></strictClockBranchRates>"
                        + "            <strictClockBranchRates><rate><parameter id=\"no.rate.opt.2\" value=\"2.0\"/></rate></strictClockBranchRates>"
                        + "        </optimalTraits>"
                        + "        <strengthOfSelectionMatrix>"
                        + "            <diagonalMatrix id=\"no.rate.selection\">"
                        + "                <parameter id=\"no.rate.alpha.diag\" value=\"1.0 1.5\" lower=\"0.0 0.0\"/>"
                        + "            </diagonalMatrix>"
                        + "        </strengthOfSelectionMatrix>",
                "",
                "");
    }

    private static String ouXmlWithTreeTimeScaledTraitRate(final String implementation) {
        return xml(implementation,
                "        <optimalTraits>"
                        + "            <strictClockBranchRates><rate><parameter id=\"tree.time.opt.1\" value=\"1.0\"/></rate></strictClockBranchRates>"
                        + "            <strictClockBranchRates><rate><parameter id=\"tree.time.opt.2\" value=\"2.0\"/></rate></strictClockBranchRates>"
                        + "        </optimalTraits>"
                        + "        <strengthOfSelectionMatrix>"
                        + "            <diagonalMatrix id=\"tree.time.selection\">"
                        + "                <parameter id=\"tree.time.alpha.diag\" value=\"1.0 1.5\" lower=\"0.0 0.0\"/>"
                        + "            </diagonalMatrix>"
                        + "        </strengthOfSelectionMatrix>",
                "    <strictClockBranchRates id=\"trait.clock.base\">"
                        + "        <rate><parameter id=\"trait.clock.base.rate\" value=\"1.0\"/></rate>"
                        + "    </strictClockBranchRates>"
                        + "    <scaledByTreeTimeBranchRates id=\"trait.clock\">"
                        + "        <treeModel idref=\"treeModel\"/>"
                        + "        <strictClockBranchRates idref=\"trait.clock.base\"/>"
                        + "        <parameter id=\"trait.clock.mean\" value=\"1.0\"/>"
                        + "    </scaledByTreeTimeBranchRates>",
                "        <scaledByTreeTimeBranchRates idref=\"trait.clock\"/>");
    }

    private static String heterogeneousOuXml(final String implementation) {
        return xml(implementation,
                "        <optimalTraits>"
                        + "            <strictClockBranchRates><rate><parameter id=\"het.opt.1\" value=\"1.0\"/></rate></strictClockBranchRates>"
                        + "            <strictClockBranchRates><rate><parameter id=\"het.opt.2\" value=\"2.0\"/></rate></strictClockBranchRates>"
                        + "        </optimalTraits>"
                        + "        <canonicalOUBranchProcesses>"
                        + "            <branchProcessIndex>"
                        + "                <parameter id=\"branch.process.index\" value=\"0 1 0 1 0\"/>"
                        + "            </branchProcessIndex>"
                        + "            <process>"
                        + "                <strengthOfSelectionMatrix>"
                        + "                    <diagonalMatrix id=\"selection.het.0\">"
                        + "                        <parameter id=\"alpha.het.0\" value=\"1.0 1.5\" lower=\"0.0 0.0\"/>"
                        + "                    </diagonalMatrix>"
                        + "                </strengthOfSelectionMatrix>"
                        + "            </process>"
                        + "            <process>"
                        + "                <strengthOfSelectionMatrix>"
                        + "                    <diagonalMatrix id=\"selection.het.1\">"
                        + "                        <parameter id=\"alpha.het.1\" value=\"2.0 0.7\" lower=\"0.0 0.0\"/>"
                        + "                    </diagonalMatrix>"
                        + "                </strengthOfSelectionMatrix>"
                        + "            </process>"
                        + "        </canonicalOUBranchProcesses>");
    }

    private static String xml(final String implementation, final String processBlock) {
        return xml(implementation,
                processBlock,
                "    <strictClockBranchRates id=\"trait.clock\">"
                        + "        <rate><parameter id=\"trait.clock.rate\" value=\"1.0\"/></rate>"
                        + "    </strictClockBranchRates>",
                "        <strictClockBranchRates idref=\"trait.clock\"/>");
    }

    private static String xml(final String implementation,
                              final String processBlock,
                              final String traitRateDefinition,
                              final String traitRateReference) {
        return "<beast>"
                + "    <taxa id=\"taxa\">"
                + "        <taxon id=\"A\"><attr name=\"X\">11 12</attr></taxon>"
                + "        <taxon id=\"B\"><attr name=\"X\">1 2</attr></taxon>"
                + "        <taxon id=\"C\"><attr name=\"X\">3 3</attr></taxon>"
                + "    </taxa>"
                + "    <newick id=\"tree\">((A:1,B:1):1,C:2);</newick>"
                + "    <treeModel id=\"treeModel\">"
                + "        <newick idref=\"tree\"/>"
                + "        <rootHeight><parameter id=\"treeModel.rootHeight\"/></rootHeight>"
                + "        <nodeHeights internalNodes=\"true\"><parameter id=\"treeModel.internalNodeHeights\"/></nodeHeights>"
                + "        <nodeHeights internalNodes=\"true\" rootNode=\"true\"><parameter id=\"treeModel.allInternalNodeHeights\"/></nodeHeights>"
                + "    </treeModel>"
                + "    <multivariateDiffusionModel id=\"diffusionModel\">"
                + "        <precisionMatrix>"
                + "            <matrixParameter id=\"precisionMatrix\">"
                + "                <parameter id=\"col1\" value=\"1.00 0.25\"/>"
                + "                <parameter id=\"col2\" value=\"0.25 1.50\"/>"
                + "            </matrixParameter>"
                + "        </precisionMatrix>"
                + "    </multivariateDiffusionModel>"
                + traitRateDefinition
                + "    <traitDataLikelihood id=\"likelihood\" traitName=\"X\" implementation=\"" + implementation
                + "\" reconstructTraits=\"false\">"
                + "        <multivariateDiffusionModel idref=\"diffusionModel\"/>"
                + "        <treeModel idref=\"treeModel\"/>"
                + traitRateReference
                + "        <traitParameter><parameter id=\"leafTraits\"/></traitParameter>"
                + "        <conjugateRootPrior>"
                + "            <meanParameter><parameter value=\"0.0 0.0\"/></meanParameter>"
                + "            <priorSampleSize><parameter value=\"10.0\"/></priorSampleSize>"
                + "        </conjugateRootPrior>"
                + processBlock
                + "    </traitDataLikelihood>"
                + "</beast>";
    }

    private static String withGradientElement(final String xml, final String gradientElement) {
        return xml.replace("</beast>", gradientElement + "</beast>");
    }
}

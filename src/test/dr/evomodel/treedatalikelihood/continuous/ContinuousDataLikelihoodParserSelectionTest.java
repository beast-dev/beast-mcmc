package test.dr.evomodel.treedatalikelihood.continuous;

import dr.app.beast.BeastParser;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.xml.XMLObject;
import junit.framework.TestCase;

import java.io.StringReader;
import java.lang.reflect.Field;
import java.util.Map;

public final class ContinuousDataLikelihoodParserSelectionTest extends TestCase {

    private static final String CDI_BM_ENGINE =
            "dr.evomodel.treedatalikelihood.continuous.CdiBMContinuousTreeLikelihoodEngine";
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
        assertFalse("BM should not select canonical OU", delegate.usesCanonicalOULikelihood());
        assertSelectedEngine(delegate, CDI_BM_ENGINE);
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

    private static ContinuousDataLikelihoodDelegate parseDelegate(final String xml) throws Exception {
        final BeastParser parser = new BeastParser(new String[0], null, false, true, true, null);
        final Map<String, XMLObject> objects = parser.parse(new StringReader(xml), true);
        final XMLObject likelihoodObject = objects.get("likelihood");
        assertNotNull("likelihood object should be stored by id", likelihoodObject);

        final TreeDataLikelihood likelihood = (TreeDataLikelihood) likelihoodObject.getNativeObject();
        return (ContinuousDataLikelihoodDelegate) likelihood.getDataLikelihoodDelegate();
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
                + "    <traitDataLikelihood id=\"likelihood\" traitName=\"X\" implementation=\"" + implementation
                + "\" reconstructTraits=\"false\">"
                + "        <multivariateDiffusionModel idref=\"diffusionModel\"/>"
                + "        <treeModel idref=\"treeModel\"/>"
                + "        <traitParameter><parameter id=\"leafTraits\"/></traitParameter>"
                + "        <conjugateRootPrior>"
                + "            <meanParameter><parameter value=\"0.0 0.0\"/></meanParameter>"
                + "            <priorSampleSize><parameter value=\"10.0\"/></priorSampleSize>"
                + "        </conjugateRootPrior>"
                + processBlock
                + "    </traitDataLikelihood>"
                + "</beast>";
    }
}

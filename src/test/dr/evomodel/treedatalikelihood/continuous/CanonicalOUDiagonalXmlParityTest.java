package test.dr.evomodel.treedatalikelihood.continuous;

import dr.app.beast.BeastParser;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.cdi.SafeMultivariateDiagonalActualizedWithDriftIntegrator;
import dr.evomodel.treedatalikelihood.continuousNew.CtsCompatibility;
import dr.xml.XMLObject;
import junit.framework.TestCase;

import java.io.File;
import java.io.FileReader;
import java.lang.reflect.Field;
import java.util.Map;

public final class CanonicalOUDiagonalXmlParityTest extends TestCase {

    private static final String XML_PATH = "tests/TestXML/testCanonicalOUDiagonalParity.xml";
    private static final String CANONICAL_OU_ENGINE =
            "dr.evomodel.treedatalikelihood.continuousNew.CanonicalOUContinuousTreeLikelihoodEngine";

    public void testXmlDiagonalOULikelihoodMatchesLegacyCDI() throws Exception {
        final ParsedLikelihoods parsed = parseLikelihoods();

        final double legacy = parsed.legacy.getLogLikelihood();
        final double canonical = parsed.canonical.getLogLikelihood();

        assertEquals("canonical-message diagonal OU XML likelihood",
                legacy,
                canonical,
                1.0e-8);
    }

    public void testXmlDiagonalOUBackendsAreSelectedAsExpected() throws Exception {
        final ParsedLikelihoods parsed = parseLikelihoods();
        final CtsCompatibility legacyDelegate =
                (CtsCompatibility) parsed.legacy.getDataLikelihoodDelegate();
        final CtsCompatibility canonicalDelegate =
                (CtsCompatibility) parsed.canonical.getDataLikelihoodDelegate();

        assertFalse("legacy XML should not select canonical OU",
                legacyDelegate.usesCanonicalOULikelihood());
        assertTrue("legacy XML should use diagonal OU CDI integrator",
                legacyDelegate.getIntegrator()
                        instanceof SafeMultivariateDiagonalActualizedWithDriftIntegrator);

        assertTrue("canonical-message XML should select canonical OU",
                canonicalDelegate.usesCanonicalOULikelihood());
        assertEquals("canonical-message XML should install canonical OU engine",
                CANONICAL_OU_ENGINE,
                privateFieldValue(canonicalDelegate, "likelihoodEngine").getClass().getName());
        assertFalse("legacy XML should not expose canonical post-order provider",
                legacyDelegate.hasCanonicalPostOrderMessageProvider());
        assertFalse("legacy XML should not expose canonical pre-order provider",
                legacyDelegate.hasCanonicalPreOrderMessageProvider());
        assertTrue("canonical-message XML should expose canonical post-order provider",
                canonicalDelegate.hasCanonicalPostOrderMessageProvider());
        assertTrue("canonical-message XML should expose canonical pre-order provider",
                canonicalDelegate.hasCanonicalPreOrderMessageProvider());
    }

    private static ParsedLikelihoods parseLikelihoods() throws Exception {
        final BeastParser parser = new BeastParser(new String[0], null, false, true, true, null);
        final Map<String, XMLObject> objects = parser.parse(new FileReader(testXmlFile()), true);

        final TreeDataLikelihood legacy = treeLikelihood(objects, "legacyDiagonalOULikelihood");
        final TreeDataLikelihood canonical = treeLikelihood(objects, "canonicalDiagonalOULikelihood");

        return new ParsedLikelihoods(legacy, canonical);
    }

    private static TreeDataLikelihood treeLikelihood(final Map<String, XMLObject> objects,
                                                     final String id) {
        final XMLObject object = objects.get(id);
        assertNotNull("XML object should be stored by id: " + id, object);
        assertTrue("XML object should be a TreeDataLikelihood: " + id,
                object.getNativeObject() instanceof TreeDataLikelihood);
        return (TreeDataLikelihood) object.getNativeObject();
    }

    private static File testXmlFile() {
        final File file = new File(XML_PATH);
        assertTrue("missing TestXML file: " + file.getAbsolutePath(), file.isFile());
        return file;
    }

    private static Object privateFieldValue(final Object target, final String fieldName) throws Exception {
        final Field field = declaredField(target.getClass(), fieldName);
        field.setAccessible(true);
        return field.get(target);
    }

    private static Field declaredField(final Class<?> initialType,
                                       final String fieldName) throws NoSuchFieldException {
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

    private static final class ParsedLikelihoods {
        final TreeDataLikelihood legacy;
        final TreeDataLikelihood canonical;

        private ParsedLikelihoods(final TreeDataLikelihood legacy,
                                  final TreeDataLikelihood canonical) {
            this.legacy = legacy;
            this.canonical = canonical;
        }
    }
}

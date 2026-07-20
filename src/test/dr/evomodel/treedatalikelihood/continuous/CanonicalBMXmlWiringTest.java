package test.dr.evomodel.treedatalikelihood.continuous;

import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.continuous.ContinuousTraitDataModel;
import dr.evomodel.treedatalikelihood.continuous.ContinuousTraitPartialsProvider;
import dr.evomodel.treedatalikelihood.continuous.DiffusionProcessDelegate;
import dr.evomodel.treedatalikelihood.continuous.HomogeneousDiffusionModelDelegate;
import dr.evomodel.treedatalikelihood.continuous.cdi.PrecisionType;
import dr.evomodel.treedatalikelihood.continuousNew.ContinuousDataLikelihoodDelegateNew;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.adapter.CanonicalBMTreeLikelihoodIntegrator;

import java.lang.reflect.Field;

public final class CanonicalBMXmlWiringTest extends ContinuousTraitTest {

    public CanonicalBMXmlWiringTest(final String name) {
        super(name);
    }

    public void testCanonicalBMLikelihoodMatchesLegacyCDI() {
        final BMSetup setup = buildSetup("bmXmlParity");

        final double legacy = setup.legacyLikelihood.getLogLikelihood();
        final double canonical = setup.canonicalLikelihood.getLogLikelihood();

        assertTrue("moment-form BM mode should be enabled", setup.canonicalDelegate.usesMomentBMLikelihood());
        assertEquals("canonical BM likelihood", legacy, canonical, 1.0e-8);
    }

    public void testCanonicalBMSelectionUsesCdiBackendOnly() {
        final BMSetup setup = buildSetup("bmXmlCdiEngineGuard");
        final Object engine = privateFieldValue(setup.canonicalDelegate, "likelihoodEngine");

        assertTrue("moment-form BM mode should be enabled", setup.canonicalDelegate.usesMomentBMLikelihood());
        assertNotNull("canonical BM selection should install a likelihood engine", engine);
        assertEquals(
                "canonical BM production path should use the native moment-form BM engine",
                "dr.evomodel.treedatalikelihood.continuousNew.MomentBMContinuousTreeLikelihoodEngine",
                engine.getClass().getName());
        assertFalse(
                "canonical BM production path should not select the natural-parameter BM integrator",
                CanonicalBMTreeLikelihoodIntegrator.class.isInstance(engine));
        assertFalse(
                "canonical BM delegate should not retain a natural-parameter BM integrator",
                hasFieldValueAssignableTo(setup.canonicalDelegate, CanonicalBMTreeLikelihoodIntegrator.class));
        assertFalse(
                "CDI-backed BM engine should not wrap a natural-parameter BM integrator",
                hasFieldValueAssignableTo(engine, CanonicalBMTreeLikelihoodIntegrator.class));
    }

    public void testCanonicalBMRefreshesTipObservationsAfterDataChangeAndRestore() {
        final BMSetup setup = buildSetup("bmXmlTipRefresh");
        final int observedTraitIndex = 0;
        final double original = traitParameter.getParameterValue(observedTraitIndex);
        final double baseline = setup.canonicalLikelihood.getLogLikelihood();

        setup.canonicalLikelihood.storeModelState();

        traitParameter.setParameterValue(observedTraitIndex, original + 0.75);
        final double changed = setup.canonicalLikelihood.getLogLikelihood();
        assertTrue("canonical BM likelihood should respond to observed tip-data changes",
                Math.abs(changed - baseline) > 1.0e-8);

        setup.canonicalLikelihood.restoreModelState();
        setup.canonicalLikelihood.makeDirty();
        final double restored = setup.canonicalLikelihood.getLogLikelihood();
        assertEquals("canonical BM likelihood should recompute from restored tip observations",
                baseline, restored, 1.0e-8);
    }

    private BMSetup buildSetup(final String tag) {
        // The native moment-form BM engine (the current "canonical BM" production backend) only
        // supports fully-observed, scalar-precision tip data (see
        // ContinuousTreeLikelihoodEngineFactory.momentBM). The inherited ContinuousTraitTest fixture
        // uses PrecisionType.FULL with real missing data (for legacy CDI factor/missing-data tests),
        // so build a fully-observed, scalar-precision view of the same underlying trait parameter
        // here and use it for both delegates, keeping the parity comparison meaningful.
        final boolean[] noMissingIndicators = new boolean[traitParameter.getDimension()];
        final ContinuousTraitPartialsProvider scalarDataModel = new ContinuousTraitDataModel(
                "scalarDataModel." + tag, traitParameter, noMissingIndicators, true,
                dimTrait, PrecisionType.SCALAR);

        final DiffusionProcessDelegate legacyProcess =
                new HomogeneousDiffusionModelDelegate(treeModel, diffusionModel);
        final ContinuousDataLikelihoodDelegate legacyDelegate =
                new ContinuousDataLikelihoodDelegate(
                        treeModel,
                        legacyProcess,
                        scalarDataModel,
                        rootPrior,
                        rateTransformation,
                        rateModel,
                        false,
                        false);
        final TreeDataLikelihood legacyLikelihood =
                new TreeDataLikelihood(legacyDelegate, treeModel, rateModel);

        final DiffusionProcessDelegate canonicalProcess =
                new HomogeneousDiffusionModelDelegate(treeModel, diffusionModel);
        final ContinuousDataLikelihoodDelegateNew canonicalDelegate =
                new ContinuousDataLikelihoodDelegateNew(
                        treeModel,
                        canonicalProcess,
                        scalarDataModel,
                        rootPrior,
                        rateTransformation,
                        rateModel,
                        false,
                        false,
                        false);
        canonicalDelegate.enableMomentBMLikelihood();
        final TreeDataLikelihood canonicalLikelihood =
                new TreeDataLikelihood(canonicalDelegate, treeModel, rateModel);

        return new BMSetup(
                tag,
                legacyDelegate,
                legacyLikelihood,
                canonicalDelegate,
                canonicalLikelihood);
    }

    private static Object privateFieldValue(final Object target, final String fieldName) {
        try {
            final Field field = declaredField(target.getClass(), fieldName);
            field.setAccessible(true);
            return field.get(target);
        } catch (IllegalAccessException e) {
            throw new AssertionError("Could not inspect field " + fieldName + " on " + target.getClass().getName(), e);
        }
    }

    private static boolean hasFieldValueAssignableTo(final Object target, final Class<?> fieldType) {
        Class<?> type = target.getClass();
        while (type != null) {
            final Field[] fields = type.getDeclaredFields();
            for (Field field : fields) {
                try {
                    field.setAccessible(true);
                    final Object value = field.get(target);
                    if (value != null && fieldType.isInstance(value)) {
                        return true;
                    }
                } catch (IllegalAccessException e) {
                    throw new AssertionError("Could not inspect field " + field.getName()
                            + " on " + target.getClass().getName(), e);
                }
            }
            type = type.getSuperclass();
        }
        return false;
    }

    private static Field declaredField(final Class<?> initialType, final String fieldName) {
        Class<?> type = initialType;
        while (type != null) {
            try {
                return type.getDeclaredField(fieldName);
            } catch (NoSuchFieldException ignored) {
                type = type.getSuperclass();
            }
        }
        throw new AssertionError("Could not find field " + fieldName + " on " + initialType.getName());
    }

    private static final class BMSetup {
        final String tag;
        final ContinuousDataLikelihoodDelegate legacyDelegate;
        final TreeDataLikelihood legacyLikelihood;
        final ContinuousDataLikelihoodDelegateNew canonicalDelegate;
        final TreeDataLikelihood canonicalLikelihood;

        private BMSetup(final String tag,
                        final ContinuousDataLikelihoodDelegate legacyDelegate,
                        final TreeDataLikelihood legacyLikelihood,
                        final ContinuousDataLikelihoodDelegateNew canonicalDelegate,
                        final TreeDataLikelihood canonicalLikelihood) {
            this.tag = tag;
            this.legacyDelegate = legacyDelegate;
            this.legacyLikelihood = legacyLikelihood;
            this.canonicalDelegate = canonicalDelegate;
            this.canonicalLikelihood = canonicalLikelihood;
        }
    }
}

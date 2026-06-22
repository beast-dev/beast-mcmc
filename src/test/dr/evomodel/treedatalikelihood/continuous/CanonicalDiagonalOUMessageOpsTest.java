package test.dr.evomodel.treedatalikelihood.continuous;

import dr.evomodel.treedatalikelihood.continuous.canonical.math.GaussianFormConverter;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianMessageOps;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianState;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianTransition;
import junit.framework.TestCase;

public final class CanonicalDiagonalOUMessageOpsTest extends TestCase {

    private static final double TOL = 1.0e-10;

    public CanonicalDiagonalOUMessageOpsTest(final String name) {
        super(name);
    }

    public void testDiagonalBackwardForwardAndTipMessagesMatchFullTransition() {
        final int d = 3;
        final double[] q = {0.82, 0.47, 0.31};
        final double[] displacement = {0.12, -0.34, 0.27};
        final double[] covariance = {
                0.90, 0.11, -0.03,
                0.11, 0.78, 0.07,
                -0.03, 0.07, 0.64
        };
        final double[] transitionMatrix = {
                q[0], 0.0, 0.0,
                0.0, q[1], 0.0,
                0.0, 0.0, q[2]
        };

        final CanonicalGaussianTransition transition = new CanonicalGaussianTransition(d);
        final GaussianFormConverter.Workspace converterWorkspace = new GaussianFormConverter.Workspace();
        converterWorkspace.ensureDim(d);
        GaussianFormConverter.fillTransitionFromMoments(
                transitionMatrix,
                displacement,
                covariance,
                d,
                converterWorkspace,
                transition);

        final CanonicalGaussianState child = new CanonicalGaussianState(d);
        GaussianFormConverter.fillStateFromMoments(
                new double[]{0.4, -0.7, 1.2},
                new double[]{
                        1.20, 0.08, 0.02,
                        0.08, 1.05, -0.04,
                        0.02, -0.04, 0.95
                },
                d,
                converterWorkspace,
                child);

        final CanonicalGaussianState parent = new CanonicalGaussianState(d);
        GaussianFormConverter.fillStateFromMoments(
                new double[]{-0.3, 0.9, 0.6},
                new double[]{
                        1.35, -0.06, 0.09,
                        -0.06, 1.10, 0.03,
                        0.09, 0.03, 1.25
                },
                d,
                converterWorkspace,
                parent);

        final CanonicalGaussianState expected = new CanonicalGaussianState(d);
        final CanonicalGaussianState actual = new CanonicalGaussianState(d);
        CanonicalGaussianMessageOps.pushBackward(
                child,
                transition,
                new CanonicalGaussianMessageOps.Workspace(d),
                expected);
        CanonicalGaussianMessageOps.pushBackwardDiagonalTransition(
                child,
                q,
                displacement,
                transition.precisionYY,
                transition.logNormalizer,
                new CanonicalGaussianMessageOps.Workspace(d),
                actual);
        assertStateClose("backward", expected, actual);

        CanonicalGaussianMessageOps.pushForward(
                parent,
                transition,
                new CanonicalGaussianMessageOps.Workspace(d),
                expected);
        CanonicalGaussianMessageOps.pushForwardDiagonalTransition(
                parent,
                q,
                displacement,
                transition.precisionYY,
                transition.logNormalizer,
                new CanonicalGaussianMessageOps.Workspace(d),
                actual);
        assertStateClose("forward", expected, actual);

        final double[] observed = {1.6, -0.2, 2.3};
        CanonicalGaussianMessageOps.conditionOnObservedSecondBlock(transition, observed, expected);
        CanonicalGaussianMessageOps.conditionOnObservedSecondBlockDiagonalTransition(
                q,
                displacement,
                transition.precisionYY,
                transition.logNormalizer,
                observed,
                actual);
        assertStateClose("exact tip", expected, actual);
    }

    private static void assertStateClose(final String label,
                                         final CanonicalGaussianState expected,
                                         final CanonicalGaussianState actual) {
        assertEquals(label + " log normalizer", expected.logNormalizer, actual.logNormalizer, TOL);
        for (int i = 0; i < expected.information.length; ++i) {
            assertEquals(label + " information[" + i + "]",
                    expected.information[i], actual.information[i], TOL);
        }
        for (int i = 0; i < expected.precision.length; ++i) {
            assertEquals(label + " precision[" + i + "]",
                    expected.precision[i], actual.precision[i], TOL);
        }
    }
}

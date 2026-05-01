package dr.evomodel.treedatalikelihood.continuous.canonical.message;

import dr.evomodel.treedatalikelihood.continuous.canonical.math.GaussianFormConverter;
import dr.evomodel.treedatalikelihood.continuous.canonical.math.MatrixOps;

/**
 * Utilities for converting expectation-form Gaussian quantities into canonical form.
 *
 * @deprecated Use {@link GaussianFormConverter} directly. That class provides the same
 *     conversions with an explicit {@link GaussianFormConverter.Workspace} parameter,
 *     avoiding the {@code ThreadLocal} allocation that this class used previously.
 *     For hot paths, pre-allocate one {@code Workspace} per long-lived computation context
 *     and reuse it across calls.
 */
@Deprecated
public final class CanonicalGaussianUtils {

    private CanonicalGaussianUtils() { }

    /**
     * @deprecated Use
     *     {@link GaussianFormConverter#fillStateFromMoments(double[], double[], int, GaussianFormConverter.Workspace, CanonicalGaussianState)}
     *     with a pre-allocated {@code Workspace}.
     */
    @Deprecated
    public static void fillStateFromMoments(final double[] mean,
                                            final double[][] covariance,
                                            final CanonicalGaussianState out) {
        final int dim = mean.length;
        final GaussianFormConverter.Workspace ws = new GaussianFormConverter.Workspace();
        ws.ensureDim(dim);
        MatrixOps.toFlat(covariance, ws.scratch, dim);
        GaussianFormConverter.fillStateFromMoments(mean, ws.scratch, dim, ws, out);
    }

    /**
     * @deprecated Use
     *     {@link GaussianFormConverter#fillTransitionFromMoments(double[], double[], double[], int, GaussianFormConverter.Workspace, CanonicalGaussianTransition)}
     *     with a pre-allocated {@code Workspace}.
     */
    @Deprecated
    public static void fillTransitionFromMoments(final double[][] transitionMatrix,
                                                 final double[] transitionOffset,
                                                 final double[][] transitionCovariance,
                                                 final CanonicalGaussianTransition out) {
        GaussianFormConverter.fillTransitionFromMoments(
                transitionMatrix, transitionOffset, transitionCovariance, out);
    }

    /**
     * @deprecated Use
     *     {@link GaussianFormConverter#fillMomentsFromState(CanonicalGaussianState, double[], double[], int, GaussianFormConverter.Workspace)}
     *     with a pre-allocated {@code Workspace}.
     */
    @Deprecated
    public static void fillMomentsFromCanonical(final CanonicalGaussianState state,
                                                final double[] meanOut,
                                                final double[][] covarianceOut) {
        GaussianFormConverter.fillMomentsFromState(state, meanOut, covarianceOut);
    }
}

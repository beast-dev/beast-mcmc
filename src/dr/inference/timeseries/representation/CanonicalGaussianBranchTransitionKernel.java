package dr.inference.timeseries.representation;

import dr.evomodel.treedatalikelihood.continuous.gaussian.CanonicalGaussianState;
import dr.evomodel.treedatalikelihood.continuous.gaussian.CanonicalGaussianTransition;

/**
 * Branch-length API for linear Gaussian transitions in canonical form.
 *
 * <p>This lives alongside {@link GaussianBranchTransitionKernel}. The existing
 * time-series and tree code still use expectation-form routines; this interface
 * provides the parallel canonical representation needed to grow a precision-form
 * pathway without disturbing the working expectation-based engines.
 */
public interface CanonicalGaussianBranchTransitionKernel {

    int getStateDimension();

    void fillInitialCanonicalState(CanonicalGaussianState out);

    void fillCanonicalTransition(double dt, CanonicalGaussianTransition out);
}

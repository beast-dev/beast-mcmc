package dr.evomodel.continuous.ou;

import dr.evomodel.treedatalikelihood.continuous.gaussian.CanonicalGaussianTransition;
import dr.inference.model.MatrixParameterInterface;

/**
 * Capability for native canonical transition construction with reusable
 * branch-local preparation.
 */
public interface CanonicalPreparedTransitionCapability {

    void fillTransitionCovariance(MatrixParameterInterface diffusionMatrix,
                                  double dt,
                                  double[][] out);

    void fillCanonicalTransition(MatrixParameterInterface diffusionMatrix,
                                 double[] stationaryMean,
                                 double dt,
                                 CanonicalGaussianTransition out);

    CanonicalPreparedBranchHandle createPreparedBranchHandle();

    CanonicalBranchWorkspace createBranchWorkspace();

    void prepareBranch(double dt,
                       double[] stationaryMean,
                       CanonicalPreparedBranchHandle prepared);

    void fillCanonicalTransitionPrepared(CanonicalPreparedBranchHandle prepared,
                                         MatrixParameterInterface diffusionMatrix,
                                         CanonicalBranchWorkspace workspace,
                                         CanonicalGaussianTransition out);
}

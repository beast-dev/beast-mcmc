package dr.evomodel.treedatalikelihood.continuous;

import dr.evolution.tree.NodeRef;
import dr.evomodel.treedatalikelihood.continuous.cdi.ContinuousDiffusionIntegrator;
import dr.evomodel.treedatalikelihood.continuous.cdi.SafeMultivariateActualizedWithDriftIntegrator;
import dr.evomodel.treedatalikelihood.preorder.BranchSufficientStatistics;
import org.ejml.data.DenseMatrix64F;

/**
 * Strategy interface for how OU actualization and its adjoint behave.
 *
 * @deprecated This interface and its implementations in {@code OUActualizationStrategies}
 *     are dead code — no external caller invokes them. The CDI path uses the separate
 *     {@link dr.evomodel.treedatalikelihood.continuous.cdi.OUActualizationStrategy}
 *     interface in the {@code cdi} package. The canonical OU path uses
 *     {@link dr.evomodel.continuous.ou.OUProcessModel} directly and does not require
 *     this strategy layer. Scheduled for removal once the legacy CDI delegate
 *     decomposition (Phase 9/10) is complete.
 */
@Deprecated
public interface OUActualizationStrategy {

    void setDiffusionStationaryVariance(SafeMultivariateActualizedWithDriftIntegrator integrator,
                                        int precisionIndex,
                                        double[] basisD,
                                        double[] basisRotations);

    void updateOrnsteinUhlenbeckDiffusionMatrices(SafeMultivariateActualizedWithDriftIntegrator integrator,
                                                  int precisionIndex,
                                                  int[] probabilityIndices,
                                                  double[] edgeLengths,
                                                  double[] optimalRates,
                                                  double[] basisD,
                                                  double[] basisRotations,
                                                  int updateCount);

    void updateIntegratedOrnsteinUhlenbeckDiffusionMatrices(SafeMultivariateActualizedWithDriftIntegrator integrator,
                                                            int precisionIndex,
                                                            int[] probabilityIndices,
                                                            double[] edgeLengths,
                                                            double[] optimalRates,
                                                            double[] basisD,
                                                            double[] basisRotations,
                                                            int updateCount);

    // Called only for non-root nodes.
    DenseMatrix64F gradientVarianceWrtVariance(NodeRef node,
                                               ContinuousDiffusionIntegrator cdi,
                                               ContinuousDataLikelihoodDelegate likelihoodDelegate,
                                               DenseMatrix64F gradient);

    // Called only for non-root nodes.
    DenseMatrix64F gradientVarianceWrtAttenuation(NodeRef node,
                                                  ContinuousDiffusionIntegrator cdi,
                                                  BranchSufficientStatistics statistics,
                                                  DenseMatrix64F dSigma);

    // Always called for internal branches (not root); vector/matrix shape matches delegate.
    DenseMatrix64F gradientDisplacementWrtDrift(NodeRef node,
                                                ContinuousDiffusionIntegrator cdi,
                                                DenseMatrix64F gradient);

    // Called only for non-root nodes.
    DenseMatrix64F gradientDisplacementWrtAttenuation(NodeRef node,
                                                      ContinuousDiffusionIntegrator cdi,
                                                      BranchSufficientStatistics statistics,
                                                      DenseMatrix64F gradient);

    // Used for heritability-style summaries.
    void meanTipVariances(double priorSampleSize,
                          double[] treeLengths,
                          DenseMatrix64F traitVariance,
                          DenseMatrix64F varSum);

    /** Root-gradient propagation: (I - A) or A depending on variant. */
    double[] rootGradient(int index,
                          ContinuousDiffusionIntegrator cdi,
                          DenseMatrix64F gradient);


    /**
     * ADDITIONS TO OUActualizationStrategy INTERFACE
     *
     * Add these method signatures to the interface.
     */

    /**
     * Compute variance gradient using direct backpropagation.
     *
     * This computes the contribution to ∂L/∂S from the variance paths:
     * - Path 1: Through A in the upward message precision Q = A^T V^{-1} A
     * - Path 3: Through V via A (stationary A-path)
     * - Path 4: Through V via Σ_stat (stationary V-path)
     *
     * @param node Current node
     * @param cdi Integrator
     * @param statistics Branch statistics
     * @param dL_dJ Gradient ∂L/∂J from root
     * @param dL_deta Gradient ∂L/∂η from root
     * @param dL_dc Gradient ∂L/∂c from root
     * @return Packed gradient (variance contribution)
     */
//    DenseMatrix64F gradientVarianceDirectBackprop(
//            NodeRef node,
//            ContinuousDiffusionIntegrator cdi,
//            BranchSufficientStatistics statistics,
//            DenseMatrix64F dL_dJ,
//            DenseMatrix64F dL_deta,
//            double dL_dc);
//
//    /**
//     * Compute displacement gradient using direct backpropagation.
//     *
//     * This computes the contribution to ∂L/∂S from the displacement path:
//     * - Path 2: Through drift b = (I - A)μ
//     *
//     * @param node Current node
//     * @param cdi Integrator
//     * @param statistics Branch statistics
//     * @param dL_dJ Gradient ∂L/∂J from root
//     * @param dL_deta Gradient ∂L/∂η from root
//     * @param dL_dc Gradient ∂L/∂c from root
//     * @return Packed gradient (displacement contribution)
//     */
//    DenseMatrix64F gradientDisplacementDirectBackprop(
//            NodeRef node,
//            ContinuousDiffusionIntegrator cdi,
//            BranchSufficientStatistics statistics,
//            DenseMatrix64F dL_dJ,
//            DenseMatrix64F dL_deta,
//            double dL_dc);
}

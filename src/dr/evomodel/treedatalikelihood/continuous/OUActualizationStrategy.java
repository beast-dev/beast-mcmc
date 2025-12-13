package dr.evomodel.treedatalikelihood.continuous;

import dr.evolution.tree.NodeRef;
import dr.evomodel.treedatalikelihood.continuous.cdi.ContinuousDiffusionIntegrator;
import dr.evomodel.treedatalikelihood.preorder.BranchSufficientStatistics;
import org.ejml.data.DenseMatrix64F;

/**
 * Strategy interface for how OU actualization and its adjoint behave.
 * The delegate is responsible for:
 *  - handling root vs non-root logic,
 *  - fixed-root logic for displacement wrt root,
 * and then forwarding to these methods.
 */
public interface OUActualizationStrategy {

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
    DenseMatrix64F gradientVarianceDirectBackprop(
            NodeRef node,
            ContinuousDiffusionIntegrator cdi,
            BranchSufficientStatistics statistics,
            DenseMatrix64F dL_dJ,
            DenseMatrix64F dL_deta,
            double dL_dc);

    /**
     * Compute displacement gradient using direct backpropagation.
     *
     * This computes the contribution to ∂L/∂S from the displacement path:
     * - Path 2: Through drift b = (I - A)μ
     *
     * @param node Current node
     * @param cdi Integrator
     * @param statistics Branch statistics
     * @param dL_dJ Gradient ∂L/∂J from root
     * @param dL_deta Gradient ∂L/∂η from root
     * @param dL_dc Gradient ∂L/∂c from root
     * @return Packed gradient (displacement contribution)
     */
    DenseMatrix64F gradientDisplacementDirectBackprop(
            NodeRef node,
            ContinuousDiffusionIntegrator cdi,
            BranchSufficientStatistics statistics,
            DenseMatrix64F dL_dJ,
            DenseMatrix64F dL_deta,
            double dL_dc);
}


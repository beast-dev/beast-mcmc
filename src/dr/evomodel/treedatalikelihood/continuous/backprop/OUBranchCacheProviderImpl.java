package dr.evomodel.treedatalikelihood.continuous.backprop;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.continuous.OUDiffusionModelDelegate;
import dr.evomodel.treedatalikelihood.continuous.cdi.ContinuousDiffusionIntegrator;
import dr.evomodel.treedatalikelihood.preorder.BranchSufficientStatistics;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

import java.util.HashMap;
import java.util.Map;

/**
 * Default implementation of OUBranchCacheProvider.
 *
 * Responsibilities:
 *  - Read branch-wise quantities (A, V) and global quantities (Σ_stat, S, μ)
 *    from the diffusion integrator and OUDiffusionModelDelegate.
 *  - Build the per-branch OU cache used by backprop gradients:
 *
 *      A      : branch actualization  (Q(t))
 *      V      : branch variance
 *      V^{-1} : inverse branch variance
 *      Σ_stat : stationary covariance (global, in original basis)
 *      b      : (I - A) μ
 *      r      : y - b, where y is the "observation" for the branch
 *      S      : selection matrix
 *      μ      : optimal value on the branch
 *      t      : branch length
 *
 *  - Cache these per node to avoid recomputation within a likelihood/gradient eval.
 */
/*
* @author Filippo Monti
 */
public class OUBranchCacheProviderImpl
        implements ContinuousTraitBackpropGradient.OUBranchCacheProvider {

    private final Tree tree;
    private final OUDiffusionModelDelegate diffusionDelegate;
    private final ContinuousDiffusionIntegrator cdi;

    private final Map<Integer, ContinuousTraitBackpropGradient.OUBranchCache> cacheMap;

    public OUBranchCacheProviderImpl(Tree tree,
                                     OUDiffusionModelDelegate diffusionDelegate,
                                     ContinuousDataLikelihoodDelegate likelihoodDelegate) {
        this.tree = tree;
        this.diffusionDelegate = diffusionDelegate;
        this.cdi = likelihoodDelegate.getIntegrator();
        this.cacheMap = new HashMap<>();
    }

    @Override
    public ContinuousTraitBackpropGradient.OUBranchCache getBranchCache(
            BranchSufficientStatistics stats, NodeRef node) {

        int nodeIndex = node.getNumber();

        ContinuousTraitBackpropGradient.OUBranchCache cache = cacheMap.get(nodeIndex);
        if (cache == null) {
            cache = computeBranchCache(stats, node);
            cacheMap.put(nodeIndex, cache);
        }

        return cache;
    }

    private ContinuousTraitBackpropGradient.OUBranchCache computeBranchCache(
            BranchSufficientStatistics stats, NodeRef node) {

        final int dim = diffusionDelegate.getDim();
        final double t = tree.getBranchLength(node);

        // ------------------------------------------------------------------
        // Global quantities
        // ------------------------------------------------------------------

        // S (selection strength) in EJML form
        DenseMatrix64F S = diffusionDelegate.getSelectionMatrixAsDense();

        // Σ_stat (stationary covariance) in original basis from CDI
        DenseMatrix64F sigmaStat = diffusionDelegate.getStationaryCovarianceFromCDI(cdi);

        // ------------------------------------------------------------------
        // Branch-specific quantities from CDI
        // ------------------------------------------------------------------

        // A = Q(t) (actualization / transition matrix)
        DenseMatrix64F A = diffusionDelegate.getBranchTransitionMatrix(node, cdi);

        // V = Var[X_child | X_parent] on this branch
        DenseMatrix64F V = diffusionDelegate.getBranchVarianceMatrix(node, cdi);

        // V^{-1}
        DenseMatrix64F Vinv = V.copy();
        CommonOps.invert(Vinv);

        // μ = optimal value on this branch (vector)
        DenseMatrix64F mu = diffusionDelegate.getOptimalValueAsDense(node);

        // ------------------------------------------------------------------
        // b = (I - A) μ
        // ------------------------------------------------------------------

        DenseMatrix64F IminusA = new DenseMatrix64F(dim, dim);
        CommonOps.setIdentity(IminusA);
        CommonOps.subtractEquals(IminusA, A);

        DenseMatrix64F b = new DenseMatrix64F(dim, 1);
        CommonOps.mult(IminusA, mu, b);

        // ------------------------------------------------------------------
        // Observation y and residual r = y - b
        // ------------------------------------------------------------------

        DenseMatrix64F y = getObservation(stats, node);

        DenseMatrix64F r = new DenseMatrix64F(dim, 1);
        CommonOps.subtract(y, b, r);

        // ------------------------------------------------------------------
        // Build cache object
        // ------------------------------------------------------------------

        return new ContinuousTraitBackpropGradient.OUBranchCache(
                A,          // Branch actualization
                V,          // Branch variance
                Vinv,       // Inverse variance
                sigmaStat,  // Stationary covariance
                b,          // (I - A) μ
                r,          // Residual y - b
                S,          // Selection matrix
                mu,         // Optimal value
                t           // Branch length
        );
    }

    /**
     * Get the "observation" for computing the residual r = y - b.
     *
     * For leaf nodes: observed tip data.
     * For internal nodes: posterior mean from aggregated child messages
     * (pseudo-observation).
     *
     * Stored in stats.getBelow().getRawMean().
     */
    private DenseMatrix64F getObservation(BranchSufficientStatistics stats, NodeRef node) {
        return stats.getBelow().getRawMean();
    }

    /** Clear all cached branch entries (e.g., when the tree / parameters change). */
    public void invalidateCache() {
        cacheMap.clear();
    }
}

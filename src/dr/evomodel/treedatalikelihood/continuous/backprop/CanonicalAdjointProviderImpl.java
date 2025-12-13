package dr.evomodel.treedatalikelihood.continuous.backprop;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.preorder.BranchConditionalDistributionDelegate;
import dr.evomodel.treedatalikelihood.preorder.BranchSufficientStatistics;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

import java.util.List;

/**
 * Canonical adjoint provider for the tree-level Gaussian root likelihood.
 *
 * In the current algebra, J_total and eta_total are simple sums of branch/root contributions,
 * hence the canonical adjoints (dL/dJ, dL/dEta, dL/dC) are the same for every node.
 *
 * This implementation therefore stores a single set of root adjoints and returns them for any nodeIndex.
 *
 * @author Filippo Monti
 */
public final class CanonicalAdjointProviderImpl
        implements ContinuousTraitBackpropGradient.CanonicalAdjointProvider {

    private static final boolean DEBUG = false;

    private final Tree tree;
    private final int dim;

    // kept for interface parity / future extensions (not used in this minimal-sums implementation)
    @SuppressWarnings("unused")
    private final ContinuousTraitBackpropGradient.OUBranchCacheProvider cacheProvider;

    private final TreeTrait<List<BranchSufficientStatistics>> treeTraitProvider;

    // Single (broadcast) canonical adjoints
    private DenseMatrix64F dL_dJ_root;     // dim x dim
    private DenseMatrix64F dL_dEta_root;   // dim x 1
    private double dL_dC_root;

    // Workspaces (allocated once, reused)
    private final DenseMatrix64F J_prior;
    private final DenseMatrix64F m_prior;
    private final DenseMatrix64F J_below;
    private final DenseMatrix64F m_below;

    private final DenseMatrix64F J_total;
    private final DenseMatrix64F eta_total;
    private final DenseMatrix64F eta_prior;

    private final DenseMatrix64F J_inv;
    private final DenseMatrix64F x_post;
    private final DenseMatrix64F x_post_outer;

    private boolean cacheValid = false;

    public CanonicalAdjointProviderImpl(String traitName,
                                        TreeDataLikelihood treeDataLikelihood,
                                        ContinuousDataLikelihoodDelegate likelihoodDelegate,
                                        ContinuousTraitBackpropGradient.OUBranchCacheProvider cacheProvider) {
        this.tree = treeDataLikelihood.getTree();
        this.dim = treeDataLikelihood.getDataLikelihoodDelegate().getTraitDim();
        this.cacheProvider = cacheProvider;

        // Ensure BranchConditionalDistributionDelegate trait exists
        String bcdName = BranchConditionalDistributionDelegate.getName(traitName);
        if (treeDataLikelihood.getTreeTrait(bcdName) == null) {
            likelihoodDelegate.addBranchConditionalDensityTrait(traitName);
        }

        @SuppressWarnings("unchecked")
        TreeTrait<List<BranchSufficientStatistics>> unchecked = treeDataLikelihood.getTreeTrait(bcdName);
        this.treeTraitProvider = unchecked;

        if (this.treeTraitProvider == null) {
            throw new IllegalStateException("TreeTrait provider for '" + bcdName + "' is null");
        }

        // Allocate workspaces once
        this.J_prior = new DenseMatrix64F(dim, dim);
        this.m_prior = new DenseMatrix64F(dim, 1);
        this.J_below = new DenseMatrix64F(dim, dim);
        this.m_below = new DenseMatrix64F(dim, 1);

        this.J_total = new DenseMatrix64F(dim, dim);
        this.eta_total = new DenseMatrix64F(dim, 1);
        this.eta_prior = new DenseMatrix64F(dim, 1);

        this.J_inv = new DenseMatrix64F(dim, dim);
        this.x_post = new DenseMatrix64F(dim, 1);
        this.x_post_outer = new DenseMatrix64F(dim, dim);
    }

    @Override
    public DenseMatrix64F get_dL_dJ(int nodeIndex) {
        ensureCacheValid();
        return dL_dJ_root;
    }

    @Override
    public DenseMatrix64F get_dL_dEta(int nodeIndex) {
        ensureCacheValid();
        return dL_dEta_root;
    }

    @Override
    public double get_dL_dC(int nodeIndex) {
        ensureCacheValid();
        return dL_dC_root;
    }

    /**
     * Compute (once) the broadcast canonical adjoints at the root:
     *  - dL/dJ_total
     *  - dL/dEta_total
     *  - dL/dC
     *
     * and store them as root-level fields that are returned for any nodeIndex.
     */
    private void ensureCacheValid() {
        if (cacheValid) return;

        final NodeRef root = tree.getRoot();

        final List<BranchSufficientStatistics> rootList = treeTraitProvider.getTrait(tree, root);
        if (rootList == null || rootList.isEmpty()) {
            throw new IllegalStateException("Root trait list is null/empty");
        }
        final BranchSufficientStatistics rootStats = rootList.get(0);

        // Grab references to upstream matrices
        final DenseMatrix64F J_prior_src = rootStats.getAbove().getRawPrecision();
        final DenseMatrix64F m_prior_src = rootStats.getAbove().getRawMean();
        final DenseMatrix64F J_below_src = rootStats.getBelow().getRawPrecision();
        final DenseMatrix64F m_below_src = rootStats.getBelow().getRawMean();

        // Copy into local workspaces (avoid aliasing/mutation risks)
        J_prior.set(J_prior_src);
        m_prior.set(m_prior_src);
        J_below.set(J_below_src);
        m_below.set(m_below_src);

        // J_total = J_below + J_prior
        J_total.set(J_below);
        CommonOps.addEquals(J_total, J_prior);

        // eta_total = J_below*m_below + J_prior*m_prior
        CommonOps.mult(J_below, m_below, eta_total);     // eta_total := eta_below
        CommonOps.mult(J_prior, m_prior, eta_prior);     // eta_prior
        CommonOps.addEquals(eta_total, eta_prior);       // eta_total += eta_prior

        if (DEBUG) {
            printMatrix("J_prior", J_prior);
            printMatrix("m_prior", m_prior);
            printMatrix("J_below", J_below);
            printMatrix("m_below", m_below);
            printMatrix("J_total", J_total);
            printMatrix("eta_total", eta_total);
        }

        // x_post = J_total^{-1} * eta_total
        J_inv.set(J_total);
        if (!CommonOps.invert(J_inv)) {
            throw new IllegalStateException("Failed to invert J_total at root");
        }
        CommonOps.mult(J_inv, eta_total, x_post);

        // dL/dJ_total = -1/2 * J^{-1} - 1/2 * x_post x_post^T
        if (dL_dJ_root == null) {
            dL_dJ_root = new DenseMatrix64F(dim, dim);
        }
        CommonOps.scale(-0.5, J_inv, dL_dJ_root);

        CommonOps.multTransB(x_post, x_post, x_post_outer);      // x_post_outer = x_post * x_post^T
        CommonOps.addEquals(dL_dJ_root, -0.5, x_post_outer);

        // dL/dEta_total = x_post
        if (dL_dEta_root == null) {
            dL_dEta_root = new DenseMatrix64F(dim, 1);
        }
        dL_dEta_root.set(x_post);

        // dL/dc = 1
        dL_dC_root = 1.0;

        if (DEBUG) {
            printMatrix("J_inv", J_inv);
            printMatrix("x_post", x_post);
            printMatrix("dL/dJ_root", dL_dJ_root);
            printMatrix("dL/dEta_root", dL_dEta_root);
            System.err.println("dL/dC_root = " + dL_dC_root);
        }

        cacheValid = true;
    }

    /**
     * Invalidate the cache when the tree/model changes.
     */
    public void invalidateCache() {
        cacheValid = false;
    }

    // --------------------------------------------------------------------
    // Debug printing
    // --------------------------------------------------------------------

    private static void printMatrix(String name, DenseMatrix64F mat) {
        if (!DEBUG) return;

        System.err.print(name);
        if (mat.numCols == 1) {
            System.err.print(" (vector): [");
            for (int i = 0; i < mat.numRows; i++) {
                System.err.printf("%.6f", mat.get(i, 0));
                if (i < mat.numRows - 1) System.err.print(", ");
            }
            System.err.println("]");
        } else {
            System.err.println(" (" + mat.numRows + "×" + mat.numCols + "):");
            for (int i = 0; i < mat.numRows; i++) {
                System.err.print("  [");
                for (int j = 0; j < mat.numCols; j++) {
                    System.err.printf("%12.6f", mat.get(i, j));
                    if (j < mat.numCols - 1) System.err.print(", ");
                }
                System.err.println("]");
            }
        }
    }
}

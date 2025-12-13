package dr.evomodel.treedatalikelihood.continuous.backprop;/*
 * ContinuousTraitBackpropGradient.java
 *
 * Backpropagation-based gradient computation for OU-type continuous
 * trait models on a single branch.
 *
 * This class coordinates the backpropagation process by:
 *   1. Fetching forward caches and canonical adjoints
 *   2. Delegating to LeafMessageBackprop for equations 28-37
 *   3. Delegating to PrimitiveParameterBackprop for PATHs 1-5
 *   4. Mapping primitive gradients to parameter space
 *
 *  BEAST is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 */

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.treedatalikelihood.continuous.ContinuousTraitGradientForBranch;
import dr.evomodel.treedatalikelihood.preorder.BranchSufficientStatistics;
import org.ejml.data.DenseMatrix64F;

/**
 * Backpropagation-based implementation of {@link ContinuousTraitGradientForBranch}.
 *
 * This class coordinates the gradient computation by delegating to specialized
 * backpropagation components:
 *   - {@link MessageBackprop}: Handles canonical → branch quantities
 *   - {@link PrimitiveParameterBackpropStrategy}: Handles branch quantities → primitive OU parameters
 */

/*
 * @author Filippo Monti
 */
public class ContinuousTraitBackpropGradient
        extends ContinuousTraitGradientForBranch.Default {

    private static final boolean DEBUG = true;

    // ------------------------------------------------------------------------
    // Dependencies
    // ------------------------------------------------------------------------

    private final OUBranchCacheProvider cacheProvider;
    private final CanonicalAdjointProvider adjointProvider;
    private final OUPrimitiveGradientMapper primitiveMapper;

    // Backpropagation components
    private final MessageBackprop leafBackprop;
    private final PrimitiveParameterBackpropStrategy primitiveStrategy;

    // ------------------------------------------------------------------------
    // Provider Interfaces
    // ------------------------------------------------------------------------

    /**
     * Provider for forward-branch OU caches (A, V, Σ_stat, b, r, etc.).
     */
    public interface OUBranchCacheProvider {
        OUBranchCache getBranchCache(BranchSufficientStatistics stats, NodeRef node);
    }



    /**
     * Provider for canonical adjoints on each branch (dℓ/dJ_i, dℓ/dη_i, dℓ/dc_i).
     */
    public interface CanonicalAdjointProvider {
        DenseMatrix64F get_dL_dJ(int nodeIndex);
        DenseMatrix64F get_dL_dEta(int nodeIndex);
        double get_dL_dC(int nodeIndex);

        default CanonicalAdjoint getAdjoint(NodeRef node) {
            int nodeIndex = node.getNumber();
            DenseMatrix64F dLdJ = get_dL_dJ(nodeIndex);
            DenseMatrix64F dLdEta = get_dL_dEta(nodeIndex);
            double dLdC = get_dL_dC(nodeIndex);
            return new CanonicalAdjoint(dLdJ, dLdEta, dLdC);
        }
    }

    /**
     * Maps primitive OU gradients to the actual parameter vector θ.
     *
     * Input:
     *   - dℓ/dS          : DenseMatrix64F (dim × dim)
     *   - dℓ/dΣ_stat     : DenseMatrix64F (dim × dim)
     *   - dℓ/dμ          : DenseMatrix64F (dim × 1)
     *
     * Output:
     *   - flattened gradient dℓ/dθ in the parameterization's coordinates.
     */
    public interface OUPrimitiveGradientMapper {
        double[] mapPrimitiveToParameters(NodeRef node,
                                          DenseMatrix64F dLdS,
                                          DenseMatrix64F dLdSigmaStat,
                                          DenseMatrix64F dLdMu,
                                          DenseMatrix64F dLdSigma);
        int getDimension();
    }

    // ------------------------------------------------------------------------
    // Data Structures
    // ------------------------------------------------------------------------

    /**
     * Forward-branch cache for OU quantities.
     * This should be filled during the forward pass by the OU integrator.
     */
    public static final class OUBranchCache {
        public final DenseMatrix64F A;
        public final DenseMatrix64F V;
        public final DenseMatrix64F Vinv;
        public final DenseMatrix64F sigmaStat;
        public final DenseMatrix64F b;
        public final DenseMatrix64F r;
        public final DenseMatrix64F S;
        public final double t;
        public final DenseMatrix64F mu;

        public OUBranchCache(DenseMatrix64F A, DenseMatrix64F V, DenseMatrix64F Vinv,
                             DenseMatrix64F sigmaStat, DenseMatrix64F b, DenseMatrix64F r,
                             DenseMatrix64F S, DenseMatrix64F mu, double t) {
            this.A = A;
            this.V = V;
            this.Vinv = Vinv;
            this.sigmaStat = sigmaStat;
            this.b = b;
            this.r = r;
            this.S = S;
            this.t = t;
            this.mu = mu;
        }
    }

    /**
     * Canonical adjoint for a branch, in the (J_i, η_i, c_i) space.
     */
    public static final class CanonicalAdjoint {
        /** dℓ/dJ_i (same shape as precision matrix) */
        public final DenseMatrix64F dLdJ;
        /** dℓ/dη_i (dim × 1) */
        public final DenseMatrix64F dLdEta;
        /** dℓ/dc_i (scalar) */
        public final double dLdC;

        public CanonicalAdjoint(DenseMatrix64F dLdJ, DenseMatrix64F dLdEta, double dLdC) {
            this.dLdJ = dLdJ;
            this.dLdEta = dLdEta;
            this.dLdC = dLdC;
        }
    }

    // ------------------------------------------------------------------------
    // Constructor
    // ------------------------------------------------------------------------

    /**
     * Create a backprop gradient computer.
     *
     * @param dim                    trait dimension
     * @param tree                   phylogenetic tree
     * @param cacheProvider          provider of forward OU branch caches
     * @param adjointProvider        provider of canonical adjoints
     * @param primitiveMapper        mapper from primitive gradients to parameters
     * @param computeSigmaGradient   if true, compute dL/dΣ in addition to dL/dΣ_stat
     */
    public ContinuousTraitBackpropGradient(int dim,
                                           Tree tree,
                                           OUBranchCacheProvider cacheProvider,
                                           CanonicalAdjointProvider adjointProvider,
                                           OUPrimitiveGradientMapper primitiveMapper,
                                           PrimitiveParameterBackpropStrategy primitiveStrategy,
                                           boolean computeSigmaGradient) {
        super(dim, tree);
        this.cacheProvider = cacheProvider;
        this.adjointProvider = adjointProvider;
        this.primitiveMapper = primitiveMapper;

        // Initialize backprop components
        this.leafBackprop = new MessageBackprop();
        this.primitiveStrategy = primitiveStrategy;
//                new PrimitiveParameterBackprop(DEBUG, computeSigmaGradient);

        if (DEBUG) {
            System.err.println("\n=== ContinuousTraitBackpropGradient Construction ===");
            System.err.println("Trait dimension: " + dim);
            System.err.println("Tree: " + tree.getId());
            System.err.println("Number of nodes: " + tree.getNodeCount());
            System.err.println("CacheProvider: " + cacheProvider.getClass().getSimpleName());
            System.err.println("AdjointProvider: " + adjointProvider.getClass().getSimpleName());
            System.err.println("PrimitiveMapper: " + primitiveMapper.getClass().getSimpleName());
            System.err.println("Compute dL/dΣ: " + computeSigmaGradient);
            System.err.println("=====================================================\n");
        }
    }

    /**
     * Convenience constructor that doesn't compute dL/dΣ (only dL/dΣ_stat).
     */
    public ContinuousTraitBackpropGradient(int dim,
                                           Tree tree,
                                           OUBranchCacheProvider cacheProvider,
                                           CanonicalAdjointProvider adjointProvider,
                                           OUPrimitiveGradientMapper primitiveMapper,
                                           PrimitiveParameterBackpropStrategy primitiveStrategy) {
        this(dim, tree, cacheProvider, adjointProvider, primitiveMapper, primitiveStrategy, false);
    }

    // ------------------------------------------------------------------------
    // Public API
    // ------------------------------------------------------------------------

    /**
     * Compute gradient w.r.t. model parameters for a branch.
     *
     * Algorithm:
     *   1. Fetch forward cache and canonical adjoint
     *   2. Backprop through leaf message (canonical → branch quantities)
     *   3. Backprop to primitive OU parameters (branch → S, Σ_stat, μ)
     *   4. Map primitive gradients to parameter vector
     *
     * @param statistics branch sufficient statistics (unused in backprop design)
     * @param node       node whose parent branch we're computing gradient for
     * @return gradient vector in parameter space
     */
    @Override
    public double[] getGradientForBranch(BranchSufficientStatistics statistics, NodeRef node) {

        // Root has no parent branch
        if (tree.isRoot(node)) {
            return new double[getDimension()];
        }

        if (DEBUG) {
            System.err.println("\n========== getGradientForBranch() ==========");
            System.err.println("Node: " + node.getNumber() + " (taxon: " +
                    (tree.isExternal(node) ? tree.getNodeTaxon(node).getId() : "internal") + ")");
            System.err.println("Branch length: " + tree.getBranchLength(node));
        }

        // Step 1: Fetch forward cache and canonical adjoints
        final OUBranchCache cache = cacheProvider.getBranchCache(statistics, node);
        final CanonicalAdjoint adjoint = adjointProvider.getAdjoint(node);

        if (cache == null || adjoint == null) {
            if (DEBUG) {
                System.err.println("WARNING: Missing cache or adjoint for node " + node.getNumber());
                System.err.println("Returning zero gradient");
                System.err.println("============================================\n");
            }
            return new double[getDimension()];
        }

        if (DEBUG) {
            printCacheContents(cache, node);
            printAdjointContents(adjoint, node);
        }

        // Step 2: Backprop through leaf message
        if (DEBUG) {
            System.err.println("\n--- Step 1: Backprop through leaf message ---");
        }
        final MessageBackprop.Result leafGrads = leafBackprop.backprop(cache, adjoint);

        if (DEBUG) {
            System.err.println("\n--- Intermediate gradients (from leaf message) ---");
            BackpropDebugPrinter.printMatrix("dL/dA", leafGrads.dLdA);
            BackpropDebugPrinter.printMatrix("dL/dV", leafGrads.dLdV);
            BackpropDebugPrinter.printMatrix("dL/dVInv", leafGrads.dLdVInv);
            BackpropDebugPrinter.printMatrix("dL/dr", leafGrads.dLdr);
            BackpropDebugPrinter.printMatrix("dL/db", leafGrads.dLdb);
        }

        // Step 3: Backprop to primitive OU parameters
        if (DEBUG) {
            System.err.println("\n--- Step 2: Backprop to primitive OU parameters ---");
        }
        PrimitiveParameterBackpropStrategy.PrimitiveGradientSet primitiveGrads =
                primitiveStrategy.backprop(cache, leafGrads);
//        final PrimitiveParameterBackprop.Result primitiveGrads =
//                primitiveBackprop.backprop(cache, leafGrads);

        if (DEBUG) {
            printPrimitiveGradients(primitiveGrads, node);
        }

        // Step 4: Map primitive gradients to parameter vector
        double[] result = primitiveMapper.mapPrimitiveToParameters(
                node,
                primitiveGrads.dLdS,
                primitiveGrads.dLdSigmaStat,
                primitiveGrads.dLdMu,
                primitiveGrads.dLdSigma
        );

        if (DEBUG) {
            System.err.println("\n--- Final Mapped Parameter Gradient ---");
            System.err.println("Dimension: " + result.length);
            BackpropDebugPrinter.printDoubleArray("Gradient", result);
            System.err.println("L2 norm: " + computeL2Norm(result));
            System.err.println("============================================\n");
        }

        return result;
    }

    @Override
    public int getParameterIndexFromNode(NodeRef node) {
        return 0;
    }

    /**
     * Dimension of parameter vector.
     *
     * TODO: This should ideally be provided by the primitiveMapper.
     * For now, hardcoded to 3 (matching typical OU parameterization).
     */
    @Override
    public int getDimension() {
        return primitiveMapper.getDimension();
    }

    // ------------------------------------------------------------------------
    // Unsupported legacy methods
    // ------------------------------------------------------------------------

    /**
     * This backprop implementation bypasses the old Q^{-1}/N chain-rule API.
     */
    @Override
    public double[] chainRule(BranchSufficientStatistics statistics,
                       NodeRef node,
                       DenseMatrix64F gradQInv,
                       DenseMatrix64F gradN) {
        throw new UnsupportedOperationException(
                "ContinuousTraitBackpropGradient does not use chainRule(gradQInv, gradN); " +
                        "call getGradientForBranch(...) directly.");
    }

    @Override
    public double[] chainRuleRoot(BranchSufficientStatistics statistics,
                           NodeRef node,
                           DenseMatrix64F gradQInv,
                           DenseMatrix64F gradN) {
        throw new UnsupportedOperationException(
                "ContinuousTraitBackpropGradient does not use chainRuleRoot(gradQInv, gradN); " +
                        "call getGradientForBranch(...) directly.");
    }

    // ------------------------------------------------------------------------
    // Debug utility methods
    // ------------------------------------------------------------------------

    private void printCacheContents(OUBranchCache cache, NodeRef node) {
        System.err.println("\n--- Forward Cache Contents (Node " + node.getNumber() + ") ---");
        System.err.println("Branch length t = " + cache.t);
        BackpropDebugPrinter.printMatrix("A = exp(-St)", cache.A);
        BackpropDebugPrinter.printMatrix("V (branch variance)", cache.V);
        BackpropDebugPrinter.printMatrix("V^{-1}", cache.Vinv);
        BackpropDebugPrinter.printMatrix("Σ_stat (stationary cov)", cache.sigmaStat);
        BackpropDebugPrinter.printMatrix("b = (I-A)μ", cache.b);
        BackpropDebugPrinter.printMatrix("r = y - b", cache.r);
        BackpropDebugPrinter.printMatrix("S (attenuation)", cache.S);
        System.err.println("---------------------------------------------------\n");
    }

    private void printAdjointContents(CanonicalAdjoint adjoint, NodeRef node) {
        System.err.println("\n--- Canonical Adjoint (Node " + node.getNumber() + ") ---");
        BackpropDebugPrinter.printMatrix("dL/dJ", adjoint.dLdJ);
        BackpropDebugPrinter.printMatrix("dL/dη", adjoint.dLdEta);
        System.err.println("dL/dc = " + adjoint.dLdC);
        System.err.println("-------------------------------------------\n");
    }

    private void printPrimitiveGradients(PrimitiveParameterBackpropStrategy.PrimitiveGradientSet pg, NodeRef node) {
        System.err.println("\n--- Primitive Gradients (Node " + node.getNumber() + ") ---");
        BackpropDebugPrinter.printMatrix("dL/dS", pg.dLdS);
        if (pg.dLdSigma != null) {
            BackpropDebugPrinter.printMatrix("dL/dΣ", pg.dLdSigma);
        }
        BackpropDebugPrinter.printMatrix("dL/dΣ_stat", pg.dLdSigmaStat);
        BackpropDebugPrinter.printMatrix("dL/dμ", pg.dLdMu);
        System.err.println("--------------------------------------------\n");
    }

    private static double computeL2Norm(double[] arr) {
        double sum = 0.0;
        for (double v : arr) {
            sum += v * v;
        }
        return Math.sqrt(sum);
    }
}
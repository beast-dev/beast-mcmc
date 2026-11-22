package dr.evomodel.treedatalikelihood.discrete;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evomodel.branchmodel.BranchModel;
import dr.evomodel.treedatalikelihood.DiscreteDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.preorder.AbstractDiscreteGradientDelegate;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.*;

/**
 * Approximate cross-product gradient for discrete substitution models.
 *
 * Computes  grad[i,j] ≈ Σ_branch  (t / L)  Σ_{c,p}  w_c w_p  q_i  p_j
 *
 * where q (pre-order) and p (post-order) are evaluated at branch bottom in
 * standard basis, L = Σ_s q_s p_s is the site likelihood, and t is the
 * effective branch length scaled by the category rate.
 *
 * Efficiency features (matching SpectralExactGradientDelegate):
 *   1. Thread-parallel branch loop via ExecutorService + per-thread workspaces.
 *   2. Single-pattern fast path: fuses the loop body for phylogeographic data
 *      (patternCount == 1), avoiding intermediate storage.
 *
 * @author Filippo Monti
 */
public final class DiscreteSubstitutionModelCrossProductDelegate extends AbstractDiscreteGradientDelegate {

    private static final String GRADIENT_TRAIT_NAME = "substitutionModelCrossProductGradient";
    private static final String THREADS_PROPERTY    = "beast.gradient.threads";

    private final String name;
    private final DiscreteDataLikelihoodDelegate likelihoodDelegate;
    private final BranchModel branchModel;
    private final int stateCount;
    private final int substitutionModelCount;

    private final int nThreads;
    private final ExecutorService pool;
    private final ThreadWorkspace[] workspaces;

    // -------------------------------------------------------------------------
    // Per-thread workspace
    // -------------------------------------------------------------------------

    private final class ThreadWorkspace {
        final double[][] accumByModel;
        final double[]   tmpPost;
        final double[]   tmpPre;

        ThreadWorkspace() {
            accumByModel = new double[substitutionModelCount][stateCount * stateCount];
            tmpPost      = new double[stateCount];
            tmpPre       = new double[stateCount];
        }
    }

    // -------------------------------------------------------------------------

    public DiscreteSubstitutionModelCrossProductDelegate(String name,
                                                         Tree tree,
                                                         DiscreteDataLikelihoodDelegate likelihoodDelegate,
                                                         int stateCount) {
        super(name, tree, likelihoodDelegate);
        this.name = name;
        this.likelihoodDelegate = likelihoodDelegate;
        this.branchModel = likelihoodDelegate.getBranchModel();
        this.stateCount = stateCount;
        this.substitutionModelCount = branchModel.getSubstitutionModels().size();

        this.nThreads = Integer.getInteger(THREADS_PROPERTY, 1);
        if (nThreads > 1) {
            this.pool = Executors.newFixedThreadPool(nThreads, r -> {
                Thread t = new Thread(r, "gradient-worker");
                t.setDaemon(true);
                return t;
            });
        } else {
            this.pool = null;
        }
        this.workspaces = new ThreadWorkspace[Math.max(nThreads, 1)];
        for (int t = 0; t < workspaces.length; t++) {
            workspaces[t] = new ThreadWorkspace();
        }
    }

    public static String getName(String name) {
        return GRADIENT_TRAIT_NAME + "." + name;
    }

    @Override
    protected int getGradientLength() {
        return stateCount * stateCount * substitutionModelCount;
    }

    @Override
    protected void getNodeDerivatives(Tree tree, double[] first, double[] second) {
        if (second != null) throw new RuntimeException("Second derivatives not yet implemented");
        if (first  == null) throw new IllegalArgumentException("First derivative buffer must not be null");
        computeGradient(first);
    }

    @Override
    protected void constructTraits(Helper treeTraitHelper) {
        treeTraitHelper.addTrait(new TreeTrait.DA() {
            @Override public String getTraitName() { return getName(name); }
            @Override public TreeTrait.Intent getIntent() { return Intent.WHOLE_TREE; }
            @Override public double[] getTrait(Tree tree, NodeRef node) { return getGradient(node); }
        });
    }

    // -------------------------------------------------------------------------

    private void computeGradient(double[] first) {
        Arrays.fill(first, 0.0);
        likelihoodDelegate.ensurePreOrderComputed();

        // Reset per-workspace accumulators
        for (ThreadWorkspace ws : workspaces) {
            for (double[] buf : ws.accumByModel) Arrays.fill(buf, 0.0);
        }

        if (nThreads <= 1) {
            processBranchRange(0, tree.getNodeCount(), workspaces[0]);
        } else {
            final int nodeCount = tree.getNodeCount();
            final int chunkSize = (nodeCount + nThreads - 1) / nThreads;
            final List<Callable<Void>> tasks = new ArrayList<>(nThreads);
            for (int t = 0; t < nThreads; t++) {
                final int start = t * chunkSize;
                final int end   = Math.min(start + chunkSize, nodeCount);
                if (start >= nodeCount) break;
                final ThreadWorkspace ws = workspaces[t];
                tasks.add(() -> { processBranchRange(start, end, ws); return null; });
            }
            try {
                for (Future<Void> f : pool.invokeAll(tasks)) f.get();
            } catch (InterruptedException | ExecutionException e) {
                throw new RuntimeException("Parallel gradient computation failed", e);
            }
        }

        // Reduce: sum per-workspace accumulators into first
        final int K2 = stateCount * stateCount;
        for (int m = 0; m < substitutionModelCount; m++) {
            final int modelOffset = m * K2;
            for (ThreadWorkspace ws : workspaces) {
                final double[] src = ws.accumByModel[m];
                for (int k = 0; k < K2; k++) {
                    first[modelOffset + k] += src[k];
                }
            }
        }
    }

    private void processBranchRange(int start, int end, ThreadWorkspace ws) {
        final double[] patternWeights  = likelihoodDelegate.getPatternWeights();
        final double[] categoryWeights = likelihoodDelegate.getCategoryWeights();
        final double[] categoryRates   = likelihoodDelegate.getSiteRateModel().getCategoryRates();

        for (int childNumber = start; childNumber < end; childNumber++) {
            final NodeRef child = tree.getNode(childNumber);
            if (tree.isRoot(child)) continue;

            final double baseBranchLength = likelihoodDelegate.getEffectiveBranchLength(childNumber);
            final BranchModel.Mapping mapping = branchModel.getBranchModelMapping(child);
            final int[]    order   = mapping.getOrder();
            final double[] weights = mapping.getWeights();

            for (int mixtureIndex = 0; mixtureIndex < order.length; mixtureIndex++) {
                final int    modelNumber = order[mixtureIndex];
                final double modelWeight = relativeWeight(mixtureIndex, weights);
                accumulateBranchContribution(
                        childNumber,
                        baseBranchLength * modelWeight,
                        categoryRates,
                        categoryWeights,
                        patternWeights,
                        ws.accumByModel[modelNumber],
                        ws.tmpPost,
                        ws.tmpPre
                );
            }
        }
    }

    private void accumulateBranchContribution(int childNumber,
                                              double branchLength,
                                              double[] categoryRates,
                                              double[] categoryWeights,
                                              double[] patternWeights,
                                              double[] destination,
                                              double[] tmpPost,
                                              double[] tmpPre) {
        final int categoryCount = likelihoodDelegate.getCategoryCount();
        final int patternCount  = likelihoodDelegate.getPatternCount();

        for (int c = 0; c < categoryCount; c++) {
            final double wc = categoryWeights[c];
            final double tc = branchLength * categoryRates[c];

            if (patternCount == 1) {
                // Fast path for single-pattern data (phylogeography): eliminates
                // the pattern loop and fuses denom computation with accumulation.
                final double wp = patternWeights[0];
                if (wp == 0.0) continue;

                likelihoodDelegate.getPostOrderBranchBottomInto(childNumber, c, 0, tmpPost);
                likelihoodDelegate.getPreOrderBranchBottomInto(childNumber, c, 0, tmpPre);

                double denom = 0.0;
                for (int s = 0; s < stateCount; s++) denom += tmpPre[s] * tmpPost[s];
                if (denom <= 0.0 || Double.isNaN(denom) || Double.isInfinite(denom)) continue;

                final double scale = (wp * wc * tc) / denom;
                for (int i = 0; i < stateCount; i++) {
                    final double qi = tmpPre[i] * scale;
                    if (qi == 0.0) continue;
                    final int rowOffset = i * stateCount;
                    for (int j = 0; j < stateCount; j++) {
                        destination[rowOffset + j] += qi * tmpPost[j];
                    }
                }
            } else {
                for (int p = 0; p < patternCount; p++) {
                    final double wp = patternWeights[p];
                    if (wp == 0.0) continue;

                    likelihoodDelegate.getPostOrderBranchBottomInto(childNumber, c, p, tmpPost);
                    likelihoodDelegate.getPreOrderBranchBottomInto(childNumber, c, p, tmpPre);

                    double denom = 0.0;
                    for (int s = 0; s < stateCount; s++) denom += tmpPre[s] * tmpPost[s];
                    if (denom <= 0.0 || Double.isNaN(denom) || Double.isInfinite(denom)) continue;

                    final double scale = (wp * wc * tc) / denom;
                    for (int i = 0; i < stateCount; i++) {
                        final double qi = tmpPre[i] * scale;
                        if (qi == 0.0) continue;
                        final int rowOffset = i * stateCount;
                        for (int j = 0; j < stateCount; j++) {
                            destination[rowOffset + j] += qi * tmpPost[j];
                        }
                    }
                }
            }
        }
    }

    private static double relativeWeight(int k, double[] weights) {
        double sum = 0.0;
        for (double w : weights) sum += w;
        return weights[k] / sum;
    }
}

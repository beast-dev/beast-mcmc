package dr.evomodel.treedatalikelihood.discrete;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evomodel.branchmodel.BranchModel;
import dr.evomodel.substmodel.EigenDecomposition;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.treedatalikelihood.DiscreteDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.discrete.beastBasedDiscreteTreeLikelihood.ComplexBlockKernelUtils;
import dr.evomodel.treedatalikelihood.discrete.beastBasedDiscreteTreeLikelihood.RealKernelUtils;
import dr.evomodel.treedatalikelihood.preorder.AbstractDiscreteGradientDelegate;

import java.util.Arrays;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.*;

/**
 * Exact spectral gradient delegate for discrete substitution models.
 *
 * Computes the exact gradient of the log-likelihood wrt the rate matrix Q
 * using Fréchet derivatives of the matrix exponential, evaluated via the
 * block-eigenvalue decomposition stored in ComplexBlockKernelUtils.
 *
 * Algorithm (from algorithmToBeImplemented.tex):
 *   For each branch i:
 *     x = R^T * preOrder_i  (rotated pre-order at parent node)
 *     y = R^{-1} * postOrder_i  (rotated post-order at child node)
 *     For each block pair (b, b'):
 *       accum += integral_0^1 e^{(1-s)t D_b^T} (x_b y_{b'}^T) e^{s t D_{b'}^T} ds
 *   gradient = R^{-T} * accum * R^T
 *
 * @author Filippo Monti
 */
public final class SpectralExactGradientDelegate extends AbstractDiscreteGradientDelegate {

    private static final String GRADIENT_TRAIT_NAME = "substitutionModelCrossProductGradient";

    /** System property to control parallelism; defaults to 1 (single-threaded). */
    private static final String THREADS_PROPERTY = "beast.gradient.threads";

    private final String name;
    private final DiscreteDataLikelihoodDelegate likelihoodDelegate;
    private final BranchModel branchModel;
    private final int stateCount;
    private final int substitutionModelCount;

    // Per-model accumulation in the rotated (eigen) basis (used only in single-thread path)
    private final double[][] eigenBasisAccumByModel;

    // Scratch for the final rotation step
    private final double[] midBuffer;

    // Parallelism
    private final int nThreads;
    private final ExecutorService pool;
    private final ThreadWorkspace[] workspaces;

    // When true, skip repeated RealKernelUtils.isAllReal() detection per gradient call.
    private final boolean forceAllReal;
    // Pre-allocated isAllReal array — reused across gradient calls to avoid hot-path allocation.
    private final boolean[] isAllRealCache;

    // -------------------------------------------------------------------------
    // Per-thread workspace
    // -------------------------------------------------------------------------

    private final class ThreadWorkspace {
        final double[][] eigenBasisAccum;
        final double[] tmpPreTop;
        final double[] tmpPreBottom;
        final double[] tmpPostBottom;
        final double[] rotatedPre;
        final double[] rotatedPost;
        final ComplexBlockKernelUtils.ComplexKernelPlan[] planByModel;
        final RealKernelUtils.RealKernelPlan[] realPlanByModel;

        ThreadWorkspace() {
            final int K2 = stateCount * stateCount;
            eigenBasisAccum  = new double[substitutionModelCount][K2];
            tmpPreTop        = new double[stateCount];
            tmpPreBottom     = new double[stateCount];
            tmpPostBottom    = new double[stateCount];
            rotatedPre       = new double[stateCount];
            rotatedPost      = new double[stateCount];
            planByModel      = new ComplexBlockKernelUtils.ComplexKernelPlan[substitutionModelCount];
            realPlanByModel  = new RealKernelUtils.RealKernelPlan[substitutionModelCount];
            for (int i = 0; i < substitutionModelCount; i++) {
                planByModel[i]     = new ComplexBlockKernelUtils.ComplexKernelPlan(stateCount);
                realPlanByModel[i] = new RealKernelUtils.RealKernelPlan(stateCount);
            }
        }
    }

    // -------------------------------------------------------------------------

    public SpectralExactGradientDelegate(String name,
                                         Tree tree,
                                         DiscreteDataLikelihoodDelegate likelihoodDelegate,
                                         int stateCount,
                                         boolean forceAllReal) {
        super(name, tree, likelihoodDelegate);
        this.name = name;
        this.likelihoodDelegate = likelihoodDelegate;
        this.branchModel = likelihoodDelegate.getBranchModel();
        this.stateCount = stateCount;
        this.substitutionModelCount = branchModel.getSubstitutionModels().size();

        final int K2 = stateCount * stateCount;
        this.eigenBasisAccumByModel = new double[substitutionModelCount][K2];
        this.midBuffer = new double[K2];

        this.forceAllReal   = forceAllReal;
        this.isAllRealCache = new boolean[substitutionModelCount];
        if (forceAllReal) {
            Arrays.fill(isAllRealCache, true);
        }

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
        if (second != null) {
            throw new RuntimeException("Second derivatives not yet implemented");
        }
        if (first == null) {
            throw new IllegalArgumentException("First derivative buffer must not be null");
        }
        computeGradient(first);
    }

    @Override
    protected void constructTraits(Helper treeTraitHelper) {
        treeTraitHelper.addTrait(new TreeTrait.DA() {
            @Override
            public String getTraitName() {
                return getName(name);
            }

            @Override
            public TreeTrait.Intent getIntent() {
                return Intent.WHOLE_TREE;
            }

            @Override
            public double[] getTrait(Tree tree, NodeRef node) {
                return getGradient(node);
            }
        });
    }

    private void computeGradient(double[] first) {
        final int K2 = stateCount * stateCount;
        Arrays.fill(first, 0.0);

        likelihoodDelegate.ensurePreOrderComputed();

        final List<SubstitutionModel> models = branchModel.getSubstitutionModels();
        final double[] patternWeights  = likelihoodDelegate.getPatternWeights();
        final double[] categoryWeights = likelihoodDelegate.getCategoryWeights();
        final double[] categoryRates   = likelihoodDelegate.getSiteRateModel().getCategoryRates();
        final boolean useInternalRotatedMessages = likelihoodDelegate.isSpectralRepresentation();

        // Cache all eigen decompositions once — avoids 598× synchronized calls per gradient eval.
        final EigenDecomposition[] eigenDecomps = new EigenDecomposition[substitutionModelCount];
        for (int m = 0; m < substitutionModelCount; m++) {
            eigenDecomps[m] = models.get(m).getEigenDecomposition();
        }

        // Detect all-real eigensystems per model (reversible CTMCs).
        // When forceAllReal=true the cache is already all-true; skip the O(K) detection calls.
        if (!forceAllReal) {
            for (int m = 0; m < substitutionModelCount; m++) {
                isAllRealCache[m] = eigenDecomps[m] != null &&
                        RealKernelUtils.isAllReal(eigenDecomps[m], stateCount);
            }
        }

        // Fill structural plan (time-independent) into every workspace.
        // Real models use RealKernelUtils (packed K² table); complex models use ComplexBlockKernelUtils.
        for (int m = 0; m < substitutionModelCount; m++) {
            if (eigenDecomps[m] != null) {
                for (ThreadWorkspace ws : workspaces) {
                    if (isAllRealCache[m]) {
                        RealKernelUtils.fillStructure(ws.realPlanByModel[m], eigenDecomps[m], stateCount);
                    } else {
                        ComplexBlockKernelUtils.fillStructure(ws.planByModel[m], eigenDecomps[m]);
                    }
                }
            }
        }

        // Reset per-workspace accumulation buffers
        for (ThreadWorkspace ws : workspaces) {
            for (double[] buf : ws.eigenBasisAccum) {
                Arrays.fill(buf, 0.0);
            }
        }

        if (nThreads <= 1) {
            processBranchRange(0, tree.getNodeCount(), workspaces[0],
                    eigenDecomps, isAllRealCache, patternWeights, categoryWeights, categoryRates, useInternalRotatedMessages);
        } else {
            final int nodeCount = tree.getNodeCount();
            final int chunkSize = (nodeCount + nThreads - 1) / nThreads;
            final List<Callable<Void>> tasks = new ArrayList<>(nThreads);
            for (int t = 0; t < nThreads; t++) {
                final int start = t * chunkSize;
                final int end   = Math.min(start + chunkSize, nodeCount);
                if (start >= nodeCount) break;
                final ThreadWorkspace ws = workspaces[t];
                tasks.add(() -> {
                    processBranchRange(start, end, ws,
                            eigenDecomps, isAllRealCache, patternWeights, categoryWeights, categoryRates, useInternalRotatedMessages);
                    return null;
                });
            }
            try {
                final List<Future<Void>> futures = pool.invokeAll(tasks);
                for (Future<Void> f : futures) f.get();
            } catch (InterruptedException | ExecutionException e) {
                throw new RuntimeException("Parallel gradient computation failed", e);
            }
        }

        // Reduce: sum per-workspace accumulations into eigenBasisAccumByModel
        for (int m = 0; m < substitutionModelCount; m++) {
            final double[] dest = eigenBasisAccumByModel[m];
            Arrays.fill(dest, 0.0);
            for (ThreadWorkspace ws : workspaces) {
                final double[] src = ws.eigenBasisAccum[m];
                for (int k = 0; k < K2; k++) {
                    dest[k] += src[k];
                }
            }
        }

        // Rotate each model's accumulation matrix and write into first
        // gradient = R^{-T} * eigenBasisAccum * R^T
        for (int m = 0; m < substitutionModelCount; m++) {
            if (eigenDecomps[m] == null) continue;
            final double[] evec = eigenDecomps[m].getEigenVectors();
            final double[] ievc = eigenDecomps[m].getInverseEigenVectors();
            rotateIntoOutput(eigenBasisAccumByModel[m], evec, ievc, first, m * K2);
        }
    }

    private void processBranchRange(int start, int end,
                                    ThreadWorkspace ws,
                                    EigenDecomposition[] eigenDecomps,
                                    boolean[] isAllReal,
                                    double[] patternWeights,
                                    double[] categoryWeights,
                                    double[] categoryRates,
                                    boolean useInternalRotatedMessages) {
        for (int childNumber = start; childNumber < end; childNumber++) {
            final NodeRef child = tree.getNode(childNumber);
            if (tree.isRoot(child)) {
                continue;
            }

            final double baseBranchLength = likelihoodDelegate.getEffectiveBranchLength(childNumber);

            final BranchModel.Mapping mapping = branchModel.getBranchModelMapping(child);
            final int[]    order   = mapping.getOrder();
            final double[] weights = mapping.getWeights();

            for (int mixtureIndex = 0; mixtureIndex < order.length; mixtureIndex++) {
                final int modelNumber    = order[mixtureIndex];
                final double modelWeight = relativeWeight(mixtureIndex, weights);
                final double branchLengthForModel = baseBranchLength * modelWeight;

                final EigenDecomposition eigenDecomp = eigenDecomps[modelNumber];
                if (eigenDecomp == null) {
                    continue;
                }

                final double[] evec = eigenDecomp.getEigenVectors();
                final double[] ievc = eigenDecomp.getInverseEigenVectors();

                if (isAllReal[modelNumber]) {
                    accumulateBranchContributionReal(
                            childNumber,
                            branchLengthForModel,
                            categoryRates,
                            categoryWeights,
                            patternWeights,
                            eigenDecomp,
                            ws.realPlanByModel[modelNumber],
                            evec,
                            ievc,
                            ws.eigenBasisAccum[modelNumber],
                            ws,
                            useInternalRotatedMessages
                    );
                } else {
                    accumulateBranchContribution(
                            childNumber,
                            branchLengthForModel,
                            categoryRates,
                            categoryWeights,
                            patternWeights,
                            eigenDecomp,
                            ws.planByModel[modelNumber],
                            evec,
                            ievc,
                            ws.eigenBasisAccum[modelNumber],
                            ws,
                            useInternalRotatedMessages
                    );
                }
            }
        }
    }

    /**
     * For a single branch (one mixture component), accumulate the Fréchet integral
     * contributions across all categories and patterns into eigenBasisAccum.
     */
    private void accumulateBranchContribution(int childNumber,
                                               double branchLength,
                                               double[] categoryRates,
                                               double[] categoryWeights,
                                               double[] patternWeights,
                                               EigenDecomposition eigenDecomp,
                                               ComplexBlockKernelUtils.ComplexKernelPlan plan,
                                               double[] evec,
                                               double[] ievc,
                                               double[] eigenBasisAccum,
                                               ThreadWorkspace ws,
                                               boolean useInternalRotatedMessages) {
        final int categoryCount = likelihoodDelegate.getCategoryCount();
        final int patternCount  = likelihoodDelegate.getPatternCount();

        for (int c = 0; c < categoryCount; c++) {
            final double wc = categoryWeights[c];
            final double tc = branchLength * categoryRates[c];

            if (patternCount == 1) {
                // Fast path for single-pattern data (e.g. phylogeography): fuse coefficient
                // computation with the outer-product accumulation in one pass, eliminating
                // the K²×16 intermediate coefficient store/load (~100 MB/gradient at K=26).
                final double wp = patternWeights[0];
                if (wp == 0.0) continue;

                double denom = 0.0;
                if (useInternalRotatedMessages) {
                    likelihoodDelegate.getInternalPreOrderBranchTopInto(childNumber, c, 0, ws.rotatedPre);
                    likelihoodDelegate.getInternalPreOrderBranchBottomInto(childNumber, c, 0, ws.tmpPreBottom);
                    likelihoodDelegate.getInternalPostOrderBranchBottomInto(childNumber, c, 0, ws.rotatedPost);
                    for (int s = 0; s < stateCount; s++) {
                        denom += ws.tmpPreBottom[s] * ws.rotatedPost[s];
                    }
                } else {
                    likelihoodDelegate.getPreOrderBranchTopInto(childNumber, c, 0, ws.tmpPreTop);
                    likelihoodDelegate.getPreOrderBranchBottomInto(childNumber, c, 0, ws.tmpPreBottom);
                    likelihoodDelegate.getPostOrderBranchBottomInto(childNumber, c, 0, ws.tmpPostBottom);
                    for (int s = 0; s < stateCount; s++) {
                        denom += ws.tmpPreBottom[s] * ws.tmpPostBottom[s];
                    }
                }
                if (denom <= 0.0 || Double.isNaN(denom) || Double.isInfinite(denom)) continue;

                if (!useInternalRotatedMessages) {
                    multiplyTransposeMatrixVector(evec, ws.tmpPreTop, ws.rotatedPre);
                    multiplyMatrixVector(ievc, ws.tmpPostBottom, ws.rotatedPost);
                }

                ComplexBlockKernelUtils.fillAndApplyToOuterProduct(
                        plan, eigenDecomp, tc, ws.rotatedPre, ws.rotatedPost,
                        (wp * wc) / denom, eigenBasisAccum);
            } else {
                // Multi-pattern path: fill coefficients once, apply across all patterns.
                ComplexBlockKernelUtils.fillTimeDependentCoefficients(plan, eigenDecomp, tc);

                for (int p = 0; p < patternCount; p++) {
                    final double wp = patternWeights[p];
                    if (wp == 0.0) continue;

                    double denom = 0.0;
                    if (useInternalRotatedMessages) {
                        likelihoodDelegate.getInternalPreOrderBranchTopInto(childNumber, c, p, ws.rotatedPre);
                        likelihoodDelegate.getInternalPreOrderBranchBottomInto(childNumber, c, p, ws.tmpPreBottom);
                        likelihoodDelegate.getInternalPostOrderBranchBottomInto(childNumber, c, p, ws.rotatedPost);
                        for (int s = 0; s < stateCount; s++) {
                            denom += ws.tmpPreBottom[s] * ws.rotatedPost[s];
                        }
                    } else {
                        likelihoodDelegate.getPreOrderBranchTopInto(childNumber, c, p, ws.tmpPreTop);
                        likelihoodDelegate.getPreOrderBranchBottomInto(childNumber, c, p, ws.tmpPreBottom);
                        likelihoodDelegate.getPostOrderBranchBottomInto(childNumber, c, p, ws.tmpPostBottom);
                        for (int s = 0; s < stateCount; s++) {
                            denom += ws.tmpPreBottom[s] * ws.tmpPostBottom[s];
                        }
                    }
                    if (denom <= 0.0 || Double.isNaN(denom) || Double.isInfinite(denom)) continue;

                    if (!useInternalRotatedMessages) {
                        multiplyTransposeMatrixVector(evec, ws.tmpPreTop, ws.rotatedPre);
                        multiplyMatrixVector(ievc, ws.tmpPostBottom, ws.rotatedPost);
                    }

                    ComplexBlockKernelUtils.applyPlanToOuterProduct(
                            plan, ws.rotatedPre, ws.rotatedPost, (wp * wc) / denom, eigenBasisAccum);
                }
            }
        }
    }

    /**
     * Same as {@link #accumulateBranchContribution} but uses the packed real-eigenvalue
     * kernel (RealKernelUtils) for models whose eigensystem is all-real.
     */
    private void accumulateBranchContributionReal(int childNumber,
                                                   double branchLength,
                                                   double[] categoryRates,
                                                   double[] categoryWeights,
                                                   double[] patternWeights,
                                                   EigenDecomposition eigenDecomp,
                                                   RealKernelUtils.RealKernelPlan plan,
                                                   double[] evec,
                                                   double[] ievc,
                                                   double[] eigenBasisAccum,
                                                   ThreadWorkspace ws,
                                                   boolean useInternalRotatedMessages) {
        final int categoryCount = likelihoodDelegate.getCategoryCount();
        final int patternCount  = likelihoodDelegate.getPatternCount();

        for (int c = 0; c < categoryCount; c++) {
            final double wc = categoryWeights[c];
            final double tc = branchLength * categoryRates[c];

            if (patternCount == 1) {
                final double wp = patternWeights[0];
                if (wp == 0.0) continue;

                double denom = 0.0;
                if (useInternalRotatedMessages) {
                    likelihoodDelegate.getInternalPreOrderBranchTopInto(childNumber, c, 0, ws.rotatedPre);
                    likelihoodDelegate.getInternalPreOrderBranchBottomInto(childNumber, c, 0, ws.tmpPreBottom);
                    likelihoodDelegate.getInternalPostOrderBranchBottomInto(childNumber, c, 0, ws.rotatedPost);
                    for (int s = 0; s < stateCount; s++) {
                        denom += ws.tmpPreBottom[s] * ws.rotatedPost[s];
                    }
                } else {
                    likelihoodDelegate.getPreOrderBranchTopInto(childNumber, c, 0, ws.tmpPreTop);
                    likelihoodDelegate.getPreOrderBranchBottomInto(childNumber, c, 0, ws.tmpPreBottom);
                    likelihoodDelegate.getPostOrderBranchBottomInto(childNumber, c, 0, ws.tmpPostBottom);
                    for (int s = 0; s < stateCount; s++) {
                        denom += ws.tmpPreBottom[s] * ws.tmpPostBottom[s];
                    }
                }
                if (denom <= 0.0 || Double.isNaN(denom) || Double.isInfinite(denom)) continue;

                if (!useInternalRotatedMessages) {
                    multiplyTransposeMatrixVector(evec, ws.tmpPreTop, ws.rotatedPre);
                    multiplyMatrixVector(ievc, ws.tmpPostBottom, ws.rotatedPost);
                }

                RealKernelUtils.fillAndApplyToOuterProduct(
                        plan, eigenDecomp, tc, ws.rotatedPre, ws.rotatedPost,
                        (wp * wc) / denom, eigenBasisAccum, stateCount);
            } else {
                RealKernelUtils.fillCoefficients(plan, eigenDecomp, tc, stateCount);

                for (int p = 0; p < patternCount; p++) {
                    final double wp = patternWeights[p];
                    if (wp == 0.0) continue;

                    double denom = 0.0;
                    if (useInternalRotatedMessages) {
                        likelihoodDelegate.getInternalPreOrderBranchTopInto(childNumber, c, p, ws.rotatedPre);
                        likelihoodDelegate.getInternalPreOrderBranchBottomInto(childNumber, c, p, ws.tmpPreBottom);
                        likelihoodDelegate.getInternalPostOrderBranchBottomInto(childNumber, c, p, ws.rotatedPost);
                        for (int s = 0; s < stateCount; s++) {
                            denom += ws.tmpPreBottom[s] * ws.rotatedPost[s];
                        }
                    } else {
                        likelihoodDelegate.getPreOrderBranchTopInto(childNumber, c, p, ws.tmpPreTop);
                        likelihoodDelegate.getPreOrderBranchBottomInto(childNumber, c, p, ws.tmpPreBottom);
                        likelihoodDelegate.getPostOrderBranchBottomInto(childNumber, c, p, ws.tmpPostBottom);
                        for (int s = 0; s < stateCount; s++) {
                            denom += ws.tmpPreBottom[s] * ws.tmpPostBottom[s];
                        }
                    }
                    if (denom <= 0.0 || Double.isNaN(denom) || Double.isInfinite(denom)) continue;

                    if (!useInternalRotatedMessages) {
                        multiplyTransposeMatrixVector(evec, ws.tmpPreTop, ws.rotatedPre);
                        multiplyMatrixVector(ievc, ws.tmpPostBottom, ws.rotatedPost);
                    }

                    RealKernelUtils.applyToOuterProduct(
                            plan, ws.rotatedPre, ws.rotatedPost, (wp * wc) / denom, eigenBasisAccum, stateCount);
                }
            }
        }
    }

    /**
     * Computes: out[modelOffset .. modelOffset+K²) = R^{-T} * eigenBasisAccum * R^T
     *
     * Step 1: mid = R^{-T} * accum
     *   mid[i,col] = sum_k R^{-T}[i,k] * accum[k,col]
     *              = sum_k ievc[k*K+i] * accum[k*K+col]
     *
     * Step 2: out = mid * R^T
     *   out[row,col] = sum_j mid[row,j] * R^T[j,col]
     *                = sum_j mid[row*K+j] * evec[col*K+j]
     */
    private void rotateIntoOutput(double[] eigenBasisAccum,
                                  double[] evec,
                                  double[] ievc,
                                  double[] out,
                                  int modelOffset) {
        final int K = stateCount;

        // Step 1: mid = R^{-T} * accum
        for (int row = 0; row < K; row++) {
            for (int col = 0; col < K; col++) {
                double sum = 0.0;
                for (int k = 0; k < K; k++) {
                    sum += ievc[k * K + row] * eigenBasisAccum[k * K + col];
                }
                midBuffer[row * K + col] = sum;
            }
        }

        // Step 2: out = mid * R^T
        for (int row = 0; row < K; row++) {
            final int rowOff = row * K;
            for (int col = 0; col < K; col++) {
                double sum = 0.0;
                final int colOff = col * K;
                for (int j = 0; j < K; j++) {
                    sum += midBuffer[rowOff + j] * evec[colOff + j];
                }
                out[modelOffset + rowOff + col] = sum;
            }
        }
    }

    private static double relativeWeight(int k, double[] weights) {
        double sum = 0.0;
        for (double w : weights) {
            sum += w;
        }
        return weights[k] / sum;
    }

    /** out[j] += matrix[i*K+j] * vector[i]  for all i  (M^T * v) */
    private void multiplyTransposeMatrixVector(double[] matrix, double[] vector, double[] out) {
        Arrays.fill(out, 0.0);
        for (int i = 0; i < stateCount; i++) {
            final double vi = vector[i];
            if (vi == 0.0) continue;
            final int rowBase = i * stateCount;
            for (int j = 0; j < stateCount; j++) {
                out[j] += matrix[rowBase + j] * vi;
            }
        }
    }

    /** out[i] = sum_j matrix[i*K+j] * vector[j]  (M * v) */
    private void multiplyMatrixVector(double[] matrix, double[] vector, double[] out) {
        for (int i = 0; i < stateCount; i++) {
            double sum = 0.0;
            final int rowBase = i * stateCount;
            for (int j = 0; j < stateCount; j++) {
                sum += matrix[rowBase + j] * vector[j];
            }
            out[i] = sum;
        }
    }
}

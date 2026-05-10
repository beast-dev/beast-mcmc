package dr.evomodel.treedatalikelihood.discrete;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evomodel.branchmodel.BranchModel;
import dr.evomodel.substmodel.EigenDecomposition;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.treedatalikelihood.DiscreteDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.discrete.beastBasedDiscreteTreeLikelihood.ComplexBlockKernelUtils;
import dr.evomodel.treedatalikelihood.preorder.AbstractDiscreteGradientDelegate;

import java.util.Arrays;
import java.util.List;

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

    private final String name;
    private final DiscreteDataLikelihoodDelegate likelihoodDelegate;
    private final BranchModel branchModel;
    private final int stateCount;
    private final int substitutionModelCount;

    // Per-branch scratch buffers
    private final double[] tmpPreTop;    // standard pre-order at parent end of branch
    private final double[] tmpPreBottom; // standard pre-order at child end (for denom)
    private final double[] tmpPostBottom; // standard post-order at child end of branch
    private final double[] rotatedPre;   // x = R^T * tmpPreTop
    private final double[] rotatedPost;  // y = R^{-1} * tmpPostBottom
    private final double[] outerProduct; // x * y^T  (stateCount x stateCount, row-major)

    // Per-model accumulation in the rotated (eigen) basis
    private final double[][] eigenBasisAccumByModel;

    // Scratch for the final rotation step
    private final double[] midBuffer;

    // Pre-allocated plan reused across all (branch, category) iterations — zero per-call allocation.
    private final ComplexBlockKernelUtils.ComplexKernelPlan planScratch;
    private final ComplexBlockKernelUtils.Workspace workspace;

    public SpectralExactGradientDelegate(String name,
                                         Tree tree,
                                         DiscreteDataLikelihoodDelegate likelihoodDelegate,
                                         int stateCount) {
        super(name, tree, likelihoodDelegate);
        this.name = name;
        this.likelihoodDelegate = likelihoodDelegate;
        this.branchModel = likelihoodDelegate.getBranchModel();
        this.stateCount = stateCount;
        this.substitutionModelCount = branchModel.getSubstitutionModels().size();

        final int K2 = stateCount * stateCount;
        this.tmpPreTop    = new double[stateCount];
        this.tmpPreBottom = new double[stateCount];
        this.tmpPostBottom = new double[stateCount];
        this.rotatedPre   = new double[stateCount];
        this.rotatedPost  = new double[stateCount];
        this.outerProduct = new double[K2];
        this.midBuffer    = new double[K2];

        this.eigenBasisAccumByModel = new double[substitutionModelCount][K2];
        this.planScratch = new ComplexBlockKernelUtils.ComplexKernelPlan(stateCount);
        this.workspace = new ComplexBlockKernelUtils.Workspace();
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

        // Reset per-model accumulation buffers
        for (double[] buf : eigenBasisAccumByModel) {
            Arrays.fill(buf, 0.0);
        }

        final List<SubstitutionModel> models = branchModel.getSubstitutionModels();
        final double[] patternWeights  = likelihoodDelegate.getPatternWeights();
        final double[] categoryWeights = likelihoodDelegate.getCategoryWeights();
        final double[] categoryRates   = likelihoodDelegate.getSiteRateModel().getCategoryRates();

        for (int childNumber = 0; childNumber < tree.getNodeCount(); childNumber++) {
            final NodeRef child = tree.getNode(childNumber);
            if (tree.isRoot(child)) {
                continue;
            }

            final double baseBranchLength = likelihoodDelegate.getEffectiveBranchLength(childNumber);

            final BranchModel.Mapping mapping = branchModel.getBranchModelMapping(child);
            final int[]    order   = mapping.getOrder();
            final double[] weights = mapping.getWeights();

            for (int mixtureIndex = 0; mixtureIndex < order.length; mixtureIndex++) {
                final int modelNumber   = order[mixtureIndex];
                final double modelWeight = relativeWeight(mixtureIndex, weights);
                final double branchLengthForModel = baseBranchLength * modelWeight;

                final EigenDecomposition eigenDecomp = models.get(modelNumber).getEigenDecomposition();
                if (eigenDecomp == null) {
                    continue;
                }

                final double[] evec = eigenDecomp.getEigenVectors();   // R
                final double[] ievc = eigenDecomp.getInverseEigenVectors(); // R^{-1}

                accumulateBranchContribution(
                        childNumber,
                        branchLengthForModel,
                        categoryRates,
                        categoryWeights,
                        patternWeights,
                        eigenDecomp,
                        evec,
                        ievc,
                        eigenBasisAccumByModel[modelNumber]
                );
            }
        }

        // Rotate each model's accumulation matrix and write into first
        // gradient = R^{-T} * eigenBasisAccum * R^T
        for (int m = 0; m < substitutionModelCount; m++) {
            final EigenDecomposition eigenDecomp = models.get(m).getEigenDecomposition();
            if (eigenDecomp == null) {
                continue;
            }
            final double[] evec = eigenDecomp.getEigenVectors();
            final double[] ievc = eigenDecomp.getInverseEigenVectors();
            final int modelOffset = m * K2;
            rotateIntoOutput(eigenBasisAccumByModel[m], evec, ievc, first, modelOffset);
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
                                               double[] evec,
                                               double[] ievc,
                                               double[] eigenBasisAccum) {
        final int categoryCount = likelihoodDelegate.getCategoryCount();
        final int patternCount  = likelihoodDelegate.getPatternCount();

        for (int c = 0; c < categoryCount; c++) {
            final double wc = categoryWeights[c];
            final double tc = branchLength * categoryRates[c];

            ComplexBlockKernelUtils.fillPlan(planScratch, eigenDecomp, tc, stateCount);
            final ComplexBlockKernelUtils.ComplexKernelPlan plan = planScratch;

            for (int p = 0; p < patternCount; p++) {
                final double wp = patternWeights[p];
                if (wp == 0.0) {
                    continue;
                }

                // Pre-order at parent end of branch (standard basis)
                likelihoodDelegate.getPreOrderBranchTopInto(childNumber, c, p, tmpPreTop);
                // Pre-order at child end (standard basis) — used only for denom normalization
                likelihoodDelegate.getPreOrderBranchBottomInto(childNumber, c, p, tmpPreBottom);
                // Post-order at child end of branch (standard basis)
                likelihoodDelegate.getPostOrderBranchBottomInto(childNumber, c, p, tmpPostBottom);

                // Denom = sum_s preBottom[s] * postBottom[s]  (= likelihood at child node)
                // Equivalent to preTop^T * P(t) * postBottom = site likelihood contribution.
                double denom = 0.0;
                for (int s = 0; s < stateCount; s++) {
                    denom += tmpPreBottom[s] * tmpPostBottom[s];
                }
                if (denom <= 0.0 || Double.isNaN(denom) || Double.isInfinite(denom)) {
                    continue;
                }

                final double scale = (wp * wc) / denom;

                // x = R^T * preTop
                multiplyTransposeMatrixVector(evec, tmpPreTop, rotatedPre);
                // y = R^{-1} * postBottom
                multiplyMatrixVector(ievc, tmpPostBottom, rotatedPost);

                // Outer product M[i][j] = scale * x[i] * y[j]
                for (int i = 0; i < stateCount; i++) {
                    final double xi = rotatedPre[i] * scale;
                    final int rowOff = i * stateCount;
                    for (int j = 0; j < stateCount; j++) {
                        outerProduct[rowOff + j] = xi * rotatedPost[j];
                    }
                }

                // Accumulate Fréchet integral contributions into eigenBasisAccum
                ComplexBlockKernelUtils.applyPlan(plan, outerProduct, eigenBasisAccum, workspace, stateCount);
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

package dr.evomodel.treedatalikelihood.discrete;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evomodel.treedatalikelihood.DiscreteDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.preorder.AbstractDiscreteGradientDelegate;

/**
 * Exact spectral gradient delegate for discrete substitution models.
 *
 * Computes the exact whole-tree gradient of the log-likelihood wrt the rate
 * matrix Q using the shared spectral branch kernel.
 *
 * @author Filippo Monti
 */
public final class SpectralExactGradientDelegate extends AbstractDiscreteGradientDelegate {

    private static final String GRADIENT_TRAIT_NAME = "substitutionModelCrossProductGradient";
    private static final String BRANCH_DIFFERENTIAL_TRAIT_NAME = "substitutionModelBranchDifferentials";
    private static final String BRANCH_LOG_RATE_SCORE_TRAIT_NAME = "substitutionModelBranchLogRateScores";

    private final String name;
    private final SpectralExactBranchKernel kernel;

    public SpectralExactGradientDelegate(String name,
                                         Tree tree,
                                         DiscreteDataLikelihoodDelegate likelihoodDelegate,
                                         int stateCount,
                                         boolean forceAllReal) {
        super(name, tree, likelihoodDelegate);
        this.name = name;
        this.kernel = new SpectralExactBranchKernel(tree, likelihoodDelegate, stateCount, forceAllReal);
    }

    public static String getName(String name) {
        return GRADIENT_TRAIT_NAME + "." + name;
    }

    public static String getBranchDifferentialTraitName(String name) {
        return BRANCH_DIFFERENTIAL_TRAIT_NAME + "." + name;
    }

    public static String getBranchLogRateScoreTraitName(String name) {
        return BRANCH_LOG_RATE_SCORE_TRAIT_NAME + "." + name;
    }

    @Override
    protected int getGradientLength() {
        return kernel.getGradientLength();
    }

    @Override
    protected void getNodeDerivatives(Tree tree, double[] first, double[] second) {
        if (second != null) {
            throw new RuntimeException("Second derivatives not yet implemented");
        }
        if (first == null) {
            throw new IllegalArgumentException("First derivative buffer must not be null");
        }
        kernel.computeGradient(first);
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
        treeTraitHelper.addTrait(new TreeTrait.DA() {
            @Override
            public String getTraitName() {
                return getBranchDifferentialTraitName(name);
            }

            @Override
            public TreeTrait.Intent getIntent() {
                return Intent.WHOLE_TREE;
            }

            @Override
            public double[] getTrait(Tree tree, NodeRef node) {
                return kernel.getBranchDifferentials();
            }
        });
        treeTraitHelper.addTrait(new TreeTrait.DA() {
            @Override
            public String getTraitName() {
                return getBranchLogRateScoreTraitName(name);
            }

            @Override
            public TreeTrait.Intent getIntent() {
                return Intent.WHOLE_TREE;
            }

            @Override
            public double[] getTrait(Tree tree, NodeRef node) {
                return kernel.getBranchLogRateScores();
            }
        });
    }
}

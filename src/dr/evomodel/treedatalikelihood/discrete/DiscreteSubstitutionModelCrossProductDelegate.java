package dr.evomodel.treedatalikelihood.discrete;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evomodel.branchmodel.BranchModel;
import dr.evomodel.treedatalikelihood.DiscreteDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.preorder.AbstractDiscreteGradientDelegate;

import java.util.Arrays;

/**
 * Discrete backend analogue of SubstitutionModelCrossProductDelegate.
 *
 * This class computes the approximate cross-product differentials directly
 * from cached pre-order and post-order messages instead of calling
 * BEAGLE.calculateCrossProductDifferentials(...).
 */
/*
 * @author Filippo Monti
 */
public final class DiscreteSubstitutionModelCrossProductDelegate extends AbstractDiscreteGradientDelegate {

    private static final String GRADIENT_TRAIT_NAME = "substitutionModelCrossProductGradient";

    private final String name;
    private final DiscreteDataLikelihoodDelegate likelihoodDelegate;
    private final BranchModel branchModel;
    private final int stateCount;
    private final int substitutionModelCount;

    private final double[] tmpPostEnd;
    private final double[] tmpPreEnd;

    public DiscreteSubstitutionModelCrossProductDelegate(String name,
                                                         Tree tree,
                                                         DiscreteDataLikelihoodDelegate likelihoodDelegate,
//                                                         BranchRateModel branchRateModel, //TODO add this too?
                                                         int stateCount) {
        super(name, tree, likelihoodDelegate);
        this.name = name;
        this.likelihoodDelegate = likelihoodDelegate;
        this.branchModel = likelihoodDelegate.getBranchModel();
        this.stateCount = stateCount;
        this.substitutionModelCount = branchModel.getSubstitutionModels().size();

        this.tmpPostEnd = new double[stateCount];
        this.tmpPreEnd = new double[stateCount];
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
        getNodeDerivatives(first);
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

    private double relativeWeight(int k, double[] weights) {
        double sum = 0.0;
        for (double w : weights) {
            sum += w;
        }
        return weights[k] / sum;
    }

    /**
     * Fills {@code first} with the approximate cross-product differentials.
     * Layout matches the BEAGLE version:
     *
     *   [ model0 KxK | model1 KxK | ... ]
     */
    private void getNodeDerivatives(double[] first) {

        final int length = stateCount * stateCount;
        Arrays.fill(first, 0.0);

        // Make sure likelihood / pre-order caches are current.
        likelihoodDelegate.ensurePreOrderComputed();

        final double[] patternWeights = likelihoodDelegate.getPatternWeights();
        final double[] categoryWeights = likelihoodDelegate.getCategoryWeights();
        final double[] categoryRates = likelihoodDelegate.getSiteRateModel().getCategoryRates();

        for (int childNumber = 0; childNumber < tree.getNodeCount(); childNumber++) {
            final NodeRef child = tree.getNode(childNumber);
            if (tree.isRoot(child)) {
                continue;
            }

            final double baseBranchLength = likelihoodDelegate.getEffectiveBranchLength(childNumber);

            final BranchModel.Mapping mapping = branchModel.getBranchModelMapping(child);
            final int[] order = mapping.getOrder();
            final double[] weights = mapping.getWeights();

            for (int mixtureIndex = 0; mixtureIndex < order.length; mixtureIndex++) {
                final int modelNumber = order[mixtureIndex];
                final double modelWeight = relativeWeight(mixtureIndex, weights);

                final int modelOffset = modelNumber * length;
                final double branchLengthForModel = baseBranchLength * modelWeight;

                accumulateBranchContribution(
                        childNumber,
                        branchLengthForModel,
                        categoryRates,
                        categoryWeights,
                        patternWeights,
                        first,
                        modelOffset
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
                                              int modelOffset) {

        final int categoryCount = likelihoodDelegate.getCategoryCount();
        final int patternCount = likelihoodDelegate.getPatternCount();

        for (int c = 0; c < categoryCount; c++) {
            final double wc = categoryWeights[c];
            final double tc = branchLength * categoryRates[c];

            for (int p = 0; p < patternCount; p++) {
                final double wp = patternWeights[p];
                if (wp == 0.0) {
                    continue;
                }

                likelihoodDelegate.getPostOrderBranchBottomInto(childNumber, c, p, tmpPostEnd);
                likelihoodDelegate.getPreOrderBranchBottomInto(childNumber, c, p, tmpPreEnd);

                double denom = 0.0;
                for (int s = 0; s < stateCount; s++) {
                    denom += tmpPreEnd[s] * tmpPostEnd[s];
                }

                if (denom <= 0.0 || Double.isNaN(denom) || Double.isInfinite(denom)) {
                    continue;
                }

                final double scale = (wp * wc * tc) / denom;

                for (int i = 0; i < stateCount; i++) {
                    final double qi = tmpPreEnd[i] * scale;
                    final int rowOffset = modelOffset + i * stateCount;
                    for (int j = 0; j < stateCount; j++) {
                        destination[rowOffset + j] += qi * tmpPostEnd[j];
                    }
                }
            }
        }
    }
}
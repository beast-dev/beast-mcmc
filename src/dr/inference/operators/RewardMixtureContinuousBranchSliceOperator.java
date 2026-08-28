package dr.inference.operators;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.branchmodel.RewardsAwareBranchModel;
import dr.evomodel.branchratemodel.RewardMixtureBranchRateModel;
import dr.evomodel.branchratemodel.RewardMixtureCategoryDecoder;
import dr.evomodel.branchratemodel.RewardMixtureCategoryDecoding;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.inference.model.Bounds;
import dr.inference.model.Parameter;
import dr.math.MathUtils;

/**
 * Branch-local slice update for total.rewards.cts on currently continuous
 * reward-mixture branches.
 *
 * This operator samples one positive-length continuous branch at a time using
 * the branch-local evidence used by the categorical Gibbs operator. It is meant
 * for reward-mixture models where this local target contains all non-constant
 * terms for the selected CTS coordinate. For dynamic reward-category ordering,
 * it also redraws the embedded category coordinate inside the continuous bucket
 * implied by the proposed CTS value.
 */
public final class RewardMixtureContinuousBranchSliceOperator
        extends SimpleMCMCOperator implements GibbsOperator {

    private static final int CONTINUOUS_CATEGORY = RewardMixtureCategoryDecoder.CONTINUOUS_CATEGORY;

    private final Parameter ctsParameter;
    private final Parameter categoryParameter;
    private final RewardMixtureCategoryDecoding categoryDecoder;
    private final RewardsAwareBranchModel rewardsAwareBranchModel;
    private final RewardsMixtureBranchWeightProvider branchWeightProvider;
    private final Tree tree;
    private final int branchCount;
    private final double windowSize;
    private final int maxSteppingOut;
    private final int maxShrinkIterations;
    private final int[] candidateBuffer;

    public RewardMixtureContinuousBranchSliceOperator(
            final Parameter categoryParameter,
            final Parameter categoryCuts,
            final RewardsAwareBranchModel rewardsAwareBranchModel,
            final TreeDataLikelihood treeDataLikelihood,
            final TreeDataLikelihood[] dependentTreeDataLikelihoods,
            final TreeDataLikelihood[] dependentContinuousTreeDataLikelihoods,
            final double windowSize,
            final int maxSteppingOut,
            final int maxShrinkIterations,
            final double weight) {
        if (categoryParameter == null) {
            throw new IllegalArgumentException("categoryParameter must be non-null");
        }
        if (categoryCuts == null) {
            throw new IllegalArgumentException("categoryCuts must be non-null");
        }
        if (rewardsAwareBranchModel == null) {
            throw new IllegalArgumentException("rewardsAwareBranchModel must be non-null");
        }
        if (!(windowSize > 0.0) || !Double.isFinite(windowSize)) {
            throw new IllegalArgumentException("windowSize must be positive and finite: " + windowSize);
        }
        if (maxSteppingOut < 1) {
            throw new IllegalArgumentException("maxSteppingOut must be at least 1: " + maxSteppingOut);
        }
        if (maxShrinkIterations < 1) {
            throw new IllegalArgumentException("maxShrinkIterations must be at least 1: " + maxShrinkIterations);
        }
        if (!(rewardsAwareBranchModel.getRateBranchModel() instanceof RewardMixtureBranchRateModel)) {
            throw new IllegalArgumentException(
                    "RewardMixtureContinuousBranchSliceOperator requires a RewardMixtureBranchRateModel");
        }

        this.rewardsAwareBranchModel = rewardsAwareBranchModel;
        this.categoryParameter = categoryParameter;
        this.ctsParameter = ((RewardMixtureBranchRateModel) rewardsAwareBranchModel.getRateBranchModel())
                .getRateParameter();
        this.branchWeightProvider = new RewardsMixtureBranchWeightProvider(
                rewardsAwareBranchModel,
                treeDataLikelihood,
                dependentTreeDataLikelihoods,
                dependentContinuousTreeDataLikelihoods);
        this.branchCount = branchWeightProvider.getBranchCount();
        if (categoryParameter.getDimension() != branchCount) {
            throw new IllegalArgumentException(
                    "categoryParameter dimension must match branch count. Found " +
                            categoryParameter.getDimension() + " but expected " + branchCount);
        }
        if (ctsParameter.getDimension() != branchCount) {
            throw new IllegalArgumentException(
                    "ctsParameter dimension must match branch count. Found " +
                            ctsParameter.getDimension() + " but expected " + branchCount);
        }

        this.categoryDecoder = resolveCategoryDecoder(rewardsAwareBranchModel, categoryParameter, categoryCuts);
        this.tree = rewardsAwareBranchModel.getTree();
        this.windowSize = windowSize;
        this.maxSteppingOut = maxSteppingOut;
        this.maxShrinkIterations = maxShrinkIterations;
        this.candidateBuffer = new int[branchCount];

        setWeight(weight);
    }

    private RewardMixtureCategoryDecoding resolveCategoryDecoder(
            final RewardsAwareBranchModel rewardsAwareBranchModel,
            final Parameter categoryParameter,
            final Parameter categoryCuts) {
        final RewardMixtureCategoryDecoding modelDecoder = rewardsAwareBranchModel.getCategoryDecoder();
        if (modelDecoder != null) {
            if (modelDecoder.getCategoryParameter() != categoryParameter ||
                    modelDecoder.getCutParameter() != categoryCuts) {
                throw new IllegalArgumentException(
                        "Slice categoryState/categoryCuts must match rewardsAwareBranchModel's category decoder");
            }
            return modelDecoder;
        }

        return new RewardMixtureCategoryDecoder(
                categoryParameter,
                categoryCuts,
                branchWeightProvider.getStateCount(),
                branchCount);
    }

    @Override
    public String getOperatorName() {
        return "RewardMixtureContinuousBranchSliceOperator(" + ctsParameter.getParameterName() + ")";
    }

    @Override
    public double doOperation() {
        refreshCategoryEmbedding();

        final int candidateCount = collectContinuousPositiveLengthBranches();
        if (candidateCount == 0) {
            return Double.NEGATIVE_INFINITY;
        }

        final int parameterIndex = candidateBuffer[MathUtils.nextInt(candidateCount)];
        branchWeightProvider.beginSingleBranchOperationCache(parameterIndex);

        final Bounds<Double> bounds = ctsParameter.getBounds();
        final double lower = bounds.getLowerLimit(parameterIndex);
        final double upper = bounds.getUpperLimit(parameterIndex);
        if (!Double.isFinite(lower) || !Double.isFinite(upper) || !(upper > lower)) {
            return Double.NEGATIVE_INFINITY;
        }

        final double current = ctsParameter.getParameterValue(parameterIndex);
        final double proposed = UnivariateSliceSampler.sample(
                current,
                lower,
                upper,
                windowSize,
                maxSteppingOut,
                maxShrinkIterations,
                new UnivariateSliceSampler.LogDensity() {
                    @Override
                    public double logDensity(final double x) {
                        return logContinuousTarget(parameterIndex, x);
                    }
                });
        if (!Double.isFinite(proposed)) {
            return Double.NEGATIVE_INFINITY;
        }

        setContinuousBranchState(parameterIndex, proposed);
        return 0.0;
    }

    private int collectContinuousPositiveLengthBranches() {
        int count = 0;
        for (int parameterIndex = 0; parameterIndex < branchCount; parameterIndex++) {
            if (categoryDecoder.getCategoryForParameterIndex(parameterIndex) != CONTINUOUS_CATEGORY) {
                continue;
            }
            final int nodeNumber = branchWeightProvider.getNodeNumberForParameterIndex(parameterIndex);
            final NodeRef node = tree.getNode(nodeNumber);
            if (tree.getBranchLength(node) > 0.0) {
                candidateBuffer[count++] = parameterIndex;
            }
        }
        return count;
    }

    private double logContinuousTarget(final int parameterIndex, final double rawReward) {
        if (!isWithinCtsBounds(parameterIndex, rawReward)) {
            return Double.NEGATIVE_INFINITY;
        }

        final double logWeight =
                branchWeightProvider.computeContinuousLogWeightForParameterIndex(parameterIndex, rawReward);
        if (!Double.isFinite(logWeight)) {
            return Double.NEGATIVE_INFINITY;
        }

        final double categoryWidth = getContinuousCategoryWidth(parameterIndex, rawReward);
        if (!(categoryWidth > 0.0) || !Double.isFinite(categoryWidth)) {
            return Double.NEGATIVE_INFINITY;
        }
        return logWeight + Math.log(categoryWidth);
    }

    private boolean isWithinCtsBounds(final int parameterIndex, final double rawReward) {
        final Bounds<Double> bounds = ctsParameter.getBounds();
        return Double.isFinite(rawReward) &&
                rawReward >= bounds.getLowerLimit(parameterIndex) &&
                rawReward <= bounds.getUpperLimit(parameterIndex);
    }

    private double getContinuousCategoryWidth(final int parameterIndex, final double rawReward) {
        return getContinuousCategoryUpperCut(parameterIndex, rawReward) -
                getContinuousCategoryLowerCut(parameterIndex, rawReward);
    }

    private double getContinuousCategoryLowerCut(final int parameterIndex, final double rawReward) {
        return categoryDecoder.getLowerCutForCategoryAtCtsValue(
                parameterIndex,
                CONTINUOUS_CATEGORY,
                rawReward);
    }

    private double getContinuousCategoryUpperCut(final int parameterIndex, final double rawReward) {
        return categoryDecoder.getUpperCutForCategoryAtCtsValue(
                parameterIndex,
                CONTINUOUS_CATEGORY,
                rawReward);
    }

    private void setContinuousBranchState(final int parameterIndex, final double rawReward) {
        final double categoryValue = sampleContinuousCategoryValue(parameterIndex, rawReward);
        ctsParameter.setParameterValueQuietly(parameterIndex, rawReward);
        categoryParameter.setParameterValueQuietly(parameterIndex, categoryValue);
        refreshCategoryEmbedding();
        ctsParameter.fireParameterChangedEvent(parameterIndex, Parameter.ChangeType.VALUE_CHANGED);
        categoryParameter.fireParameterChangedEvent(parameterIndex, Parameter.ChangeType.VALUE_CHANGED);
    }

    private double sampleContinuousCategoryValue(final int parameterIndex, final double rawReward) {
        final double lower = getContinuousCategoryLowerCut(parameterIndex, rawReward);
        final double upper = getContinuousCategoryUpperCut(parameterIndex, rawReward);
        final double width = upper - lower;
        if (!(width > 0.0) || !Double.isFinite(width)) {
            throw new IllegalStateException("Invalid continuous-category interval for branch " + parameterIndex);
        }
        double u = MathUtils.nextDouble();
        if (u == 0.0) {
            u = Math.nextUp(0.0);
        }
        return lower + u * width;
    }

    private void refreshCategoryEmbedding() {
        if (rewardsAwareBranchModel.getCategoryDecoder() == null) {
            categoryDecoder.refreshEmbedding();
        } else {
            rewardsAwareBranchModel.refreshCategoryDecoderEmbedding();
        }
    }
}

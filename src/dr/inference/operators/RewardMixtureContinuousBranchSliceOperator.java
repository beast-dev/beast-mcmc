package dr.inference.operators;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.branchmodel.RewardsAwareBranchModel;
import dr.evomodel.branchratemodel.PerBranchRewardMixtureCategoryDecoder;
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
        categoryDecoder.refreshEmbedding();
        rewardsAwareBranchModel.refreshCategoryDecoderEmbedding();

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
        final double currentLogDensity = logContinuousTarget(parameterIndex, current);
        if (!Double.isFinite(currentLogDensity)) {
            return Double.NEGATIVE_INFINITY;
        }

        final double cutoff = currentLogDensity + MathUtils.randomLogDouble();
        final SliceInterval interval = constructInterval(parameterIndex, current, cutoff, lower, upper);
        final double proposed = drawFromInterval(parameterIndex, current, cutoff, interval);
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

    private SliceInterval constructInterval(final int parameterIndex,
                                            final double current,
                                            final double cutoff,
                                            final double lower,
                                            final double upper) {
        double left = current - windowSize * MathUtils.nextDouble();
        double right = left + windowSize;
        if (left < lower) {
            left = lower;
        }
        if (right > upper) {
            right = upper;
        }

        int leftSteps = MathUtils.nextInt(maxSteppingOut);
        int rightSteps = (maxSteppingOut - 1) - leftSteps;

        while (leftSteps > 0 && left > lower) {
            final double nextLeft = Math.max(lower, left - windowSize);
            if (!(logContinuousTarget(parameterIndex, nextLeft) > cutoff)) {
                break;
            }
            left = nextLeft;
            leftSteps--;
        }

        while (rightSteps > 0 && right < upper) {
            final double nextRight = Math.min(upper, right + windowSize);
            if (!(logContinuousTarget(parameterIndex, nextRight) > cutoff)) {
                break;
            }
            right = nextRight;
            rightSteps--;
        }

        return new SliceInterval(left, right);
    }

    private double drawFromInterval(final int parameterIndex,
                                    final double current,
                                    final double cutoff,
                                    final SliceInterval interval) {
        double left = interval.left;
        double right = interval.right;
        for (int i = 0; i < maxShrinkIterations; i++) {
            final double proposed = MathUtils.uniform(left, right);
            if (logContinuousTarget(parameterIndex, proposed) >= cutoff) {
                return proposed;
            }
            if (proposed < current) {
                left = proposed;
            } else {
                right = proposed;
            }
        }
        return Double.NaN;
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
        if (categoryDecoder instanceof PerBranchRewardMixtureCategoryDecoder) {
            return ((PerBranchRewardMixtureCategoryDecoder) categoryDecoder)
                    .getLowerCutForCategoryAtCtsValue(parameterIndex, CONTINUOUS_CATEGORY, rawReward);
        }
        return categoryDecoder.getLowerCut(parameterIndex, CONTINUOUS_CATEGORY);
    }

    private double getContinuousCategoryUpperCut(final int parameterIndex, final double rawReward) {
        if (categoryDecoder instanceof PerBranchRewardMixtureCategoryDecoder) {
            return ((PerBranchRewardMixtureCategoryDecoder) categoryDecoder)
                    .getUpperCutForCategoryAtCtsValue(parameterIndex, CONTINUOUS_CATEGORY, rawReward);
        }
        return categoryDecoder.getUpperCut(parameterIndex, CONTINUOUS_CATEGORY);
    }

    private void setContinuousBranchState(final int parameterIndex, final double rawReward) {
        final double categoryValue = sampleContinuousCategoryValue(parameterIndex, rawReward);
        ctsParameter.setParameterValueQuietly(parameterIndex, rawReward);
        categoryParameter.setParameterValueQuietly(parameterIndex, categoryValue);
        categoryDecoder.refreshEmbedding();
        rewardsAwareBranchModel.refreshCategoryDecoderEmbedding();
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

    private static final class SliceInterval {
        private final double left;
        private final double right;

        private SliceInterval(final double left, final double right) {
            this.left = left;
            this.right = right;
        }
    }
}

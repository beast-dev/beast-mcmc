package dr.inference.operators;

import dr.evomodel.branchmodel.RewardMixtureAtomicPseudoPrior;
import dr.evomodel.branchmodel.RewardsAwareBranchModel;
import dr.evomodel.branchratemodel.RewardMixtureBranchRateModel;
import dr.evomodel.branchratemodel.RewardMixtureCategoryDecoder;
import dr.evomodel.branchratemodel.RewardMixtureCategoryDecoding;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.inference.model.Bounds;
import dr.inference.model.Parameter;
import dr.math.MathUtils;

import java.util.Arrays;

/**
 * Branch-local Gibbs-within-slice update for reward-mixture category and CTS
 * reward coordinates.
 *
 * The BEAST state stores the reward category as an embedded continuous
 * coordinate. This operator samples a decoded category conditional on the
 * current CTS value, samples the CTS value conditional on that decoded
 * category, and finally redraws the embedded coordinate uniformly inside the
 * decoded category interval implied by the proposed CTS value.
 */
public final class RewardMixtureBranchJointCtsCategoryOperator
        extends SimpleMCMCOperator implements GibbsOperator {

    private final Parameter ctsParameter;
    private final Parameter categoryParameter;
    private final RewardMixtureCategoryDecoding categoryDecoder;
    private final RewardsAwareBranchModel rewardsAwareBranchModel;
    private final RewardsMixtureBranchWeightProvider branchWeightProvider;
    private final RewardMixtureAtomicPseudoPrior atomicPseudoPrior;
    private final int branchCount;
    private final int categoryCount;
    private final double windowSize;
    private final int maxSteppingOut;
    private final int maxShrinkIterations;
    private final double[] logCategoryMasses;

    public RewardMixtureBranchJointCtsCategoryOperator(
            final Parameter categoryParameter,
            final Parameter categoryCuts,
            final RewardsAwareBranchModel rewardsAwareBranchModel,
            final TreeDataLikelihood treeDataLikelihood,
            final TreeDataLikelihood[] dependentTreeDataLikelihoods,
            final TreeDataLikelihood[] dependentContinuousTreeDataLikelihoods,
            final RewardMixtureAtomicPseudoPrior atomicPseudoPrior,
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
        if (!(rewardsAwareBranchModel.getRateBranchModel() instanceof RewardMixtureBranchRateModel)) {
            throw new IllegalArgumentException(
                    "RewardMixtureBranchJointCtsCategoryOperator requires a RewardMixtureBranchRateModel");
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

        this.rewardsAwareBranchModel = rewardsAwareBranchModel;
        this.categoryParameter = categoryParameter;
        this.ctsParameter = ((RewardMixtureBranchRateModel) rewardsAwareBranchModel.getRateBranchModel())
                .getRateParameter();
        this.atomicPseudoPrior = atomicPseudoPrior;
        if (atomicPseudoPrior != null && atomicPseudoPrior.getParameter() != ctsParameter) {
            throw new IllegalArgumentException(
                    "atomicPseudoPrior must target the CTS reward parameter used by rewardsAwareBranchModel");
        }

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
        this.categoryCount = categoryDecoder.getCategoryCount();
        this.windowSize = windowSize;
        this.maxSteppingOut = maxSteppingOut;
        this.maxShrinkIterations = maxShrinkIterations;
        this.logCategoryMasses = new double[categoryCount];

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
                        "Joint category/CTS categoryState/categoryCuts must match rewardsAwareBranchModel's category decoder");
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
        return "RewardMixtureBranchJointCtsCategoryOperator(" +
                ctsParameter.getParameterName() + "," +
                categoryParameter.getParameterName() + ")";
    }

    @Override
    public double doOperation() {
        refreshCategoryEmbedding();

        final int parameterIndex = MathUtils.nextInt(branchCount);
        branchWeightProvider.beginSingleBranchOperationCache(parameterIndex);

        final Bounds<Double> bounds = ctsParameter.getBounds();
        final double lower = bounds.getLowerLimit(parameterIndex);
        final double upper = bounds.getUpperLimit(parameterIndex);
        if (!Double.isFinite(lower) || !Double.isFinite(upper) || !(upper > lower)) {
            return Double.NEGATIVE_INFINITY;
        }

        final double oldCts = ctsParameter.getParameterValue(parameterIndex);
        final double oldCategory = categoryParameter.getParameterValue(parameterIndex);
        final int proposedCategory = sampleCategory(parameterIndex, oldCts);
        if (proposedCategory < 0) {
            return Double.NEGATIVE_INFINITY;
        }

        final double proposedCts = UnivariateSliceSampler.sample(
                oldCts,
                lower,
                upper,
                windowSize,
                maxSteppingOut,
                maxShrinkIterations,
                new UnivariateSliceSampler.LogDensity() {
                    @Override
                    public double logDensity(final double x) {
                        return logCtsTarget(parameterIndex, proposedCategory, x);
                    }
                });
        if (!Double.isFinite(proposedCts)) {
            restoreBranchState(parameterIndex, oldCts, oldCategory);
            return Double.NEGATIVE_INFINITY;
        }

        final double proposedEmbeddedCategory =
                sampleValueInCategoryAtCtsValue(parameterIndex, proposedCategory, proposedCts);
        if (!Double.isFinite(proposedEmbeddedCategory)) {
            restoreBranchState(parameterIndex, oldCts, oldCategory);
            return Double.NEGATIVE_INFINITY;
        }

        setBranchState(parameterIndex, proposedCts, proposedEmbeddedCategory);
        Arrays.fill(logCategoryMasses, Double.NEGATIVE_INFINITY);
        return 0.0;
    }

    private int sampleCategory(final int parameterIndex, final double rawReward) {
        final RewardsMixtureBranchResamplingHelper.BranchWeights weights =
                branchWeightProvider.getOperationCachedBranchWeightsForParameterIndex(parameterIndex);

        double logTotal = Double.NEGATIVE_INFINITY;
        for (int category = 0; category < categoryCount; category++) {
            final double logWeight = branchWeightProvider.getLogWeightForCategory(weights, category);
            final double logMass = logWeight +
                    getPseudoPriorLogDensity(parameterIndex, category, rawReward) +
                    getLogCategoryWidth(parameterIndex, category, rawReward);
            logCategoryMasses[category] = Double.isFinite(logMass)
                    ? logMass
                    : Double.NEGATIVE_INFINITY;
            logTotal = RewardsMixtureBranchResamplingHelper.logAdd(logTotal, logCategoryMasses[category]);
        }

        if (!Double.isFinite(logTotal)) {
            return -1;
        }

        final double u = MathUtils.nextDouble();
        double cumulative = 0.0;
        int lastFiniteCategory = -1;
        for (int category = 0; category < categoryCount; category++) {
            if (Double.isFinite(logCategoryMasses[category])) {
                lastFiniteCategory = category;
                cumulative += Math.exp(logCategoryMasses[category] - logTotal);
            }
            if (u < cumulative) {
                return category;
            }
        }
        return lastFiniteCategory;
    }

    private double logCtsTarget(final int parameterIndex,
                                final int category,
                                final double rawReward) {
        if (!isWithinCtsBounds(parameterIndex, rawReward)) {
            return Double.NEGATIVE_INFINITY;
        }

        final double logWeight =
                branchWeightProvider.computeLogWeightForCategoryAtValue(parameterIndex, category, rawReward);
        if (!Double.isFinite(logWeight)) {
            return Double.NEGATIVE_INFINITY;
        }

        final double logPseudoPrior = getPseudoPriorLogDensity(parameterIndex, category, rawReward);
        if (!Double.isFinite(logPseudoPrior)) {
            return Double.NEGATIVE_INFINITY;
        }

        final double logWidth = getLogCategoryWidth(parameterIndex, category, rawReward);
        if (!Double.isFinite(logWidth)) {
            return Double.NEGATIVE_INFINITY;
        }

        return logWeight + logPseudoPrior + logWidth;
    }

    private double getPseudoPriorLogDensity(final int parameterIndex,
                                            final int category,
                                            final double rawReward) {
        return atomicPseudoPrior == null
                ? 0.0
                : atomicPseudoPrior.getLogDensityForCategoryAtValue(parameterIndex, category, rawReward);
    }

    private double getLogCategoryWidth(final int parameterIndex,
                                       final int category,
                                       final double rawReward) {
        final double lower = categoryDecoder.getLowerCutForCategoryAtCtsValue(
                parameterIndex,
                category,
                rawReward);
        final double upper = categoryDecoder.getUpperCutForCategoryAtCtsValue(
                parameterIndex,
                category,
                rawReward);
        final double width = upper - lower;
        return width > 0.0 && Double.isFinite(width)
                ? Math.log(width)
                : Double.NEGATIVE_INFINITY;
    }

    private boolean isWithinCtsBounds(final int parameterIndex, final double rawReward) {
        final Bounds<Double> bounds = ctsParameter.getBounds();
        return Double.isFinite(rawReward) &&
                rawReward >= bounds.getLowerLimit(parameterIndex) &&
                rawReward <= bounds.getUpperLimit(parameterIndex);
    }

    private double sampleValueInCategoryAtCtsValue(final int parameterIndex,
                                                   final int category,
                                                   final double rawReward) {
        final double lower = categoryDecoder.getLowerCutForCategoryAtCtsValue(
                parameterIndex,
                category,
                rawReward);
        final double upper = categoryDecoder.getUpperCutForCategoryAtCtsValue(
                parameterIndex,
                category,
                rawReward);
        final double width = upper - lower;
        if (!(width > 0.0) || !Double.isFinite(width)) {
            return Double.NaN;
        }

        double u = MathUtils.nextDouble();
        if (category > RewardMixtureCategoryDecoder.CONTINUOUS_CATEGORY && u == 0.0) {
            u = Math.nextUp(0.0);
        }
        return lower + u * width;
    }

    private void setBranchState(final int parameterIndex,
                                final double rawReward,
                                final double embeddedCategoryValue) {
        ctsParameter.setParameterValueQuietly(parameterIndex, rawReward);
        categoryParameter.setParameterValueQuietly(parameterIndex, embeddedCategoryValue);
        refreshCategoryEmbedding();
        ctsParameter.fireParameterChangedEvent(parameterIndex, Parameter.ChangeType.VALUE_CHANGED);
        categoryParameter.fireParameterChangedEvent(parameterIndex, Parameter.ChangeType.VALUE_CHANGED);
    }

    private void restoreBranchState(final int parameterIndex,
                                    final double rawReward,
                                    final double embeddedCategoryValue) {
        ctsParameter.setParameterValueQuietly(parameterIndex, rawReward);
        categoryParameter.setParameterValueQuietly(parameterIndex, embeddedCategoryValue);
        refreshCategoryEmbedding();
        ctsParameter.fireParameterChangedEvent(parameterIndex, Parameter.ChangeType.VALUE_CHANGED);
        categoryParameter.fireParameterChangedEvent(parameterIndex, Parameter.ChangeType.VALUE_CHANGED);
    }

    private void refreshCategoryEmbedding() {
        if (rewardsAwareBranchModel.getCategoryDecoder() == null) {
            categoryDecoder.refreshEmbedding();
        } else {
            rewardsAwareBranchModel.refreshCategoryDecoderEmbedding();
        }
    }
}

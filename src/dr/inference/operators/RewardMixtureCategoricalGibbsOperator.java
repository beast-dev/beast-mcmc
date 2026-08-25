package dr.inference.operators;

import dr.evomodel.branchmodel.RewardsAwareBranchModel;
import dr.evomodel.branchmodel.RewardMixtureAtomicPseudoPrior;
import dr.evomodel.branchratemodel.RewardMixtureCategoryDecoder;
import dr.evomodel.branchratemodel.RewardMixtureCategoryDecoding;
import dr.evomodel.branchratemodel.RewardMixtureBranchRateModel;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.inference.model.Parameter;
import dr.math.MathUtils;

import java.util.Arrays;

/**
 * Branch-local Gibbs refresh for the embedded categorical reward-mixture state.
 *
 * The target parameter is continuous, but its value is interpreted through fixed
 * cuts: category 0 is the continuous reward branch, and categories 1..K are
 * atomic no-jump states.  Conditional on a selected decoded category, this
 * operator samples uniformly within that interval so the embedded coordinate
 * remains a genuine continuous parameter.
 *
 * @author Filippo Monti
 */
public final class RewardMixtureCategoricalGibbsOperator extends SimpleMCMCOperator implements GibbsOperator {

    private final Parameter categoryParameter;
    private final RewardMixtureCategoryDecoding categoryDecoder;
    private final RewardsMixtureBranchWeightProvider branchWeightProvider;
    private final RewardMixtureAtomicPseudoPrior atomicPseudoPrior;
    private final TreeDataLikelihood treeDataLikelihood;
    private final TreeDataLikelihood[] dependentTreeDataLikelihoods;
    private final TreeDataLikelihood[] dependentContinuousTreeDataLikelihoods;
    private final int branchCount;
    private final int categoryCount;
    private final double updateProportion;
    private final int[] candidateBuffer;
    private final double[] storedCategoryValues;
    private final double[] logCategoryMasses;

    public RewardMixtureCategoricalGibbsOperator(
            final Parameter categoryParameter,
            final Parameter categoryCuts,
            final RewardsAwareBranchModel rewardsAwareBranchModel,
            final TreeDataLikelihood treeDataLikelihood,
            final TreeDataLikelihood[] dependentTreeDataLikelihoods,
            final TreeDataLikelihood[] dependentContinuousTreeDataLikelihoods,
            final double updateProportion,
            final double weight) {
        this(categoryParameter,
                categoryCuts,
                rewardsAwareBranchModel,
                treeDataLikelihood,
                dependentTreeDataLikelihoods,
                dependentContinuousTreeDataLikelihoods,
                null,
                updateProportion,
                weight);
    }

    public RewardMixtureCategoricalGibbsOperator(
            final Parameter categoryParameter,
            final Parameter categoryCuts,
            final RewardsAwareBranchModel rewardsAwareBranchModel,
            final TreeDataLikelihood treeDataLikelihood,
            final TreeDataLikelihood[] dependentTreeDataLikelihoods,
            final TreeDataLikelihood[] dependentContinuousTreeDataLikelihoods,
            final RewardMixtureAtomicPseudoPrior atomicPseudoPrior,
            final double updateProportion,
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
        if (updateProportion <= 0.0 || updateProportion > 1.0) {
            throw new IllegalArgumentException("updateProportion must be in (0, 1]. Found: " + updateProportion);
        }
        if (atomicPseudoPrior != null &&
                rewardsAwareBranchModel.getRateBranchModel() instanceof RewardMixtureBranchRateModel) {
            final RewardMixtureBranchRateModel rewardMixtureBranchRateModel =
                    (RewardMixtureBranchRateModel) rewardsAwareBranchModel.getRateBranchModel();
            if (atomicPseudoPrior.getParameter() != rewardMixtureBranchRateModel.getRateParameter()) {
                throw new IllegalArgumentException(
                        "atomicPseudoPrior must target the cts reward parameter used by rewardsAwareBranchModel");
            }
        }

        this.categoryParameter = categoryParameter;
        this.treeDataLikelihood = treeDataLikelihood;
        this.dependentTreeDataLikelihoods = dependentTreeDataLikelihoods == null
                ? new TreeDataLikelihood[0]
                : Arrays.copyOf(dependentTreeDataLikelihoods, dependentTreeDataLikelihoods.length);
        this.dependentContinuousTreeDataLikelihoods = dependentContinuousTreeDataLikelihoods == null
                ? new TreeDataLikelihood[0]
                : Arrays.copyOf(dependentContinuousTreeDataLikelihoods, dependentContinuousTreeDataLikelihoods.length);
        this.atomicPseudoPrior = atomicPseudoPrior;
        this.branchWeightProvider = new RewardsMixtureBranchWeightProvider(
                rewardsAwareBranchModel,
                treeDataLikelihood,
                this.dependentTreeDataLikelihoods,
                this.dependentContinuousTreeDataLikelihoods);
        this.branchCount = branchWeightProvider.getBranchCount();
        if (categoryParameter.getDimension() != branchCount) {
            throw new IllegalArgumentException(
                    "categoryParameter dimension must match branch count. Found " +
                            categoryParameter.getDimension() + " but expected " + branchCount);
        }

        this.categoryDecoder = resolveCategoryDecoder(rewardsAwareBranchModel, categoryParameter, categoryCuts);
        this.categoryCount = categoryDecoder.getCategoryCount();
        this.updateProportion = updateProportion;
        this.candidateBuffer = new int[branchCount];
        this.storedCategoryValues = new double[branchCount];
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
                        "Gibbs categoryState/categoryCuts must match rewardsAwareBranchModel's category decoder");
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
        return "RewardMixtureCategoricalGibbsOperator(" + categoryParameter.getParameterName() + ")";
    }

    @Override
    public double doOperation() {
        categoryDecoder.refreshEmbedding();
        storeCurrentCategoryValues();

        final int updateCount = Math.max(1,
                Math.min(branchCount, (int) Math.round(updateProportion * branchCount)));
        initializeCandidateBuffer();
        shufflePrefix(updateCount);

        branchWeightProvider.beginOperationCache();
        boolean valid = true;
        boolean changedAnyCategory = false;

        for (int i = 0; i < updateCount; i++) {
            final int parameterIndex = candidateBuffer[i];
            final int resampleStatus = resampleBranch(parameterIndex);
            valid &= resampleStatus >= 0;

            final boolean changedCategory = resampleStatus > 0;
            changedAnyCategory |= changedCategory;
            if (valid) {
                RewardMixturePerformanceStats.recordCategoricalGibbsBranchUpdate(changedCategory);
            }
            categoryParameter.fireParameterChangedEvent(parameterIndex, Parameter.ChangeType.VALUE_CHANGED);

            if (!valid) {
                break;
            }
            if (changedCategory && i + 1 < updateCount) {
                branchWeightProvider.clearOperationCache(
                        RewardMixturePerformanceStats.OperationCacheClearReason.GIBBS_CATEGORY_CHANGE);
            } else if (!changedCategory && i + 1 < updateCount) {
                RewardMixturePerformanceStats.recordCategoricalGibbsSkippedCacheClearAfterSameCategory();
            }
        }

        if (!valid || (changedAnyCategory && !currentLogTargetIsFinite())) {
            restoreStoredCategoryValuesQuietly();
            categoryParameter.fireParameterChangedEvent();
            return Double.NEGATIVE_INFINITY;
        }

        return 0.0;
    }

    private boolean currentLogTargetIsFinite() {
        treeDataLikelihood.makeDirty();
        if (!Double.isFinite(treeDataLikelihood.getLogLikelihood())) {
            return false;
        }
        for (TreeDataLikelihood likelihood : dependentTreeDataLikelihoods) {
            likelihood.makeDirty();
            if (!Double.isFinite(likelihood.getLogLikelihood())) {
                return false;
            }
        }
        for (TreeDataLikelihood likelihood : dependentContinuousTreeDataLikelihoods) {
            likelihood.makeDirty();
            if (!Double.isFinite(likelihood.getLogLikelihood())) {
                return false;
            }
        }
        if (atomicPseudoPrior != null) {
            atomicPseudoPrior.makeDirty();
            if (!Double.isFinite(atomicPseudoPrior.getLogLikelihood())) {
                return false;
            }
        }
        return true;
    }

    /**
     * @return -1 for invalid weights, 0 for a within-category redraw, 1 for a
     * decoded-category change.
     */
    private int resampleBranch(final int parameterIndex) {
        final RewardsMixtureBranchResamplingHelper.BranchWeights weights =
                branchWeightProvider.getOperationCachedBranchWeightsForParameterIndex(parameterIndex);

        final int currentCategory = categoryDecoder.getCategoryForParameterIndex(parameterIndex);
        final int category = sampleCategory(parameterIndex, weights);
        if (category < 0) {
            return -1;
        }

        categoryParameter.setParameterValueQuietly(parameterIndex, sampleValueInCategory(parameterIndex, category));
        return category == currentCategory ? 0 : 1;
    }

    private int sampleCategory(final int parameterIndex,
                               final RewardsMixtureBranchResamplingHelper.BranchWeights weights) {
        double logTotal = Double.NEGATIVE_INFINITY;
        for (int category = 0; category < categoryCount; category++) {
            final double lower = categoryDecoder.getLowerCut(parameterIndex, category);
            final double upper = categoryDecoder.getUpperCut(parameterIndex, category);
            final double width = upper - lower;
            final double logWeight = branchWeightProvider.getLogWeightForCategory(weights, category) +
                    getPseudoPriorLogDensity(parameterIndex, category);

            if (!(width > 0.0) || !Double.isFinite(logWeight)) {
                logCategoryMasses[category] = Double.NEGATIVE_INFINITY;
            } else {
                logCategoryMasses[category] = logWeight + Math.log(width);
            }
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

    private double sampleValueInCategory(final int parameterIndex, final int category) {
        final double lower = categoryDecoder.getLowerCut(parameterIndex, category);
        final double upper = categoryDecoder.getUpperCut(parameterIndex, category);
        final double width = upper - lower;

        double u = MathUtils.nextDouble();
        if (category > 0 && u == 0.0) {
            u = Math.nextUp(0.0);
        }
        return lower + u * width;
    }

    private double getPseudoPriorLogDensity(final int parameterIndex, final int category) {
        return atomicPseudoPrior == null ? 0.0 : atomicPseudoPrior.getLogDensityForCategory(parameterIndex, category);
    }

    private void initializeCandidateBuffer() {
        for (int i = 0; i < branchCount; i++) {
            candidateBuffer[i] = i;
        }
    }

    private void shufflePrefix(final int updateCount) {
        for (int i = 0; i < updateCount; i++) {
            final int j = i + MathUtils.nextInt(branchCount - i);
            final int tmp = candidateBuffer[i];
            candidateBuffer[i] = candidateBuffer[j];
            candidateBuffer[j] = tmp;
        }
    }

    private void storeCurrentCategoryValues() {
        for (int i = 0; i < branchCount; i++) {
            storedCategoryValues[i] = categoryParameter.getParameterValue(i);
        }
    }

    private void restoreStoredCategoryValuesQuietly() {
        for (int i = 0; i < branchCount; i++) {
            categoryParameter.setParameterValueQuietly(i, storedCategoryValues[i]);
        }
        Arrays.fill(logCategoryMasses, Double.NEGATIVE_INFINITY);
    }
}

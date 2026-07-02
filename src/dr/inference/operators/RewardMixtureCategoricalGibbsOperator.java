package dr.inference.operators;

import dr.evomodel.branchmodel.RewardsAwareBranchModel;
import dr.evomodel.branchratemodel.RewardMixtureCategoryDecoder;
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
    private final RewardMixtureCategoryDecoder categoryDecoder;
    private final RewardsMixtureBranchWeightProvider branchWeightProvider;
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

        this.categoryParameter = categoryParameter;
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

        this.categoryDecoder = new RewardMixtureCategoryDecoder(
                categoryParameter,
                categoryCuts,
                branchWeightProvider.getStateCount(),
                branchCount);
        this.categoryCount = categoryDecoder.getCategoryCount();
        this.updateProportion = updateProportion;
        this.candidateBuffer = new int[branchCount];
        this.storedCategoryValues = new double[branchCount];
        this.logCategoryMasses = new double[categoryCount];

        setWeight(weight);
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

        branchWeightProvider.refreshLikelihoodMessages();
        boolean valid = true;

        for (int i = 0; i < updateCount; i++) {
            final int parameterIndex = candidateBuffer[i];
            valid &= resampleBranch(parameterIndex);
            categoryParameter.fireParameterChangedEvent(parameterIndex, Parameter.ChangeType.VALUE_CHANGED);

            if (!valid) {
                break;
            }
            if (i + 1 < updateCount) {
                branchWeightProvider.refreshLikelihoodMessages();
            }
        }

        if (!valid) {
            restoreStoredCategoryValuesQuietly();
            categoryParameter.fireParameterChangedEvent();
            return Double.NEGATIVE_INFINITY;
        }

        return 0.0;
    }

    private boolean resampleBranch(final int parameterIndex) {
        final RewardsMixtureBranchResamplingHelper.BranchWeights weights =
                branchWeightProvider.computeBranchWeightsForParameterIndex(parameterIndex);

        final int category = sampleCategory(weights);
        if (category < 0) {
            return false;
        }

        categoryParameter.setParameterValueQuietly(parameterIndex, sampleValueInCategory(category));
        return true;
    }

    private int sampleCategory(final RewardsMixtureBranchResamplingHelper.BranchWeights weights) {
        double logTotal = Double.NEGATIVE_INFINITY;
        for (int category = 0; category < categoryCount; category++) {
            final double lower = categoryDecoder.getLowerCut(category);
            final double upper = categoryDecoder.getUpperCut(category);
            final double width = upper - lower;
            final double logWeight = branchWeightProvider.getLogWeightForCategory(weights, category);

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
        for (int category = 0; category < categoryCount; category++) {
            if (Double.isFinite(logCategoryMasses[category])) {
                cumulative += Math.exp(logCategoryMasses[category] - logTotal);
            }
            if (u < cumulative) {
                return category;
            }
        }

        return categoryCount - 1;
    }

    private double sampleValueInCategory(final int category) {
        final double lower = categoryDecoder.getLowerCut(category);
        final double upper = categoryDecoder.getUpperCut(category);
        final double width = upper - lower;

        double u = MathUtils.nextDouble();
        if (category > 0 && u == 0.0) {
            u = Math.nextUp(0.0);
        }
        return lower + u * width;
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

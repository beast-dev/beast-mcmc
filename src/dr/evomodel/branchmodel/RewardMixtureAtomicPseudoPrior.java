package dr.evomodel.branchmodel;

import dr.evomodel.branchratemodel.ArbitraryBranchRates;
import dr.evomodel.branchratemodel.RewardMixtureCategoryDecoder;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.hmc.NumericGradientStepSizeProvider;
import dr.inference.model.Bounds;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.math.distributions.NormalDistribution;
import dr.xml.Reportable;

import java.util.Arrays;

/**
 * Conditional prior for the latent continuous reward coordinate in a
 * reward-mixture branch model.
 *
 * Continuous branches receive a uniform density on the cts parameter bounds.
 * Atomic branches receive a truncated-normal pseudo-prior centered on the
 * selected atomic reward. This keeps inactive cts coordinates near useful
 * values without changing the branch's atomic CTMC likelihood.
 *
 * @author Filippo Monti
 */
public final class RewardMixtureAtomicPseudoPrior extends Likelihood.Abstract
        implements GradientWrtParameterProvider, Reportable {

    public static final String MODEL_NAME = "RewardMixtureAtomicPseudoPrior";

    private static final double LOG_SQRT_TWO_PI = 0.5 * Math.log(2.0 * Math.PI);

    private final RewardsAwareBranchModel rewardsAwareBranchModel;
    private final ArbitraryBranchRates totalRewardsBranchRates;
    private final Parameter parameter;
    private final double standardDeviation;
    private final double variance;
    private final double[] gradient;

    private double numericGradientStepSize =
            NumericGradientStepSizeProvider.StepSizeLevel.MEDIUM.getStepSizeRatio();

    public RewardMixtureAtomicPseudoPrior(final RewardsAwareBranchModel rewardsAwareBranchModel,
                                          final ArbitraryBranchRates totalRewardsBranchRates,
                                          final double standardDeviation) {
        super(rewardsAwareBranchModel);

        if (rewardsAwareBranchModel == null) {
            throw new IllegalArgumentException("rewardsAwareBranchModel must be non-null");
        }
        if (totalRewardsBranchRates == null) {
            throw new IllegalArgumentException("totalRewardsBranchRates must be non-null");
        }
        if (rewardsAwareBranchModel.getRateBranchModel() != totalRewardsBranchRates) {
            throw new IllegalArgumentException(
                    "totalRewardsBranchRates must be the branch-rate model used by rewardsAwareBranchModel");
        }
        if (rewardsAwareBranchModel.getCategoryDecoder() == null) {
            throw new IllegalArgumentException(
                    "RewardMixtureAtomicPseudoPrior requires categorical reward-mixture branch states");
        }
        if (!(standardDeviation > 0.0) || !Double.isFinite(standardDeviation)) {
            throw new IllegalArgumentException("standardDeviation must be positive and finite: " + standardDeviation);
        }

        this.rewardsAwareBranchModel = rewardsAwareBranchModel;
        this.totalRewardsBranchRates = totalRewardsBranchRates;
        this.parameter = totalRewardsBranchRates.getRateParameter();
        this.standardDeviation = standardDeviation;
        this.variance = standardDeviation * standardDeviation;
        this.gradient = new double[parameter.getDimension()];

        validateDimension();
        validateBounds();
    }

    @Override
    public Likelihood getLikelihood() {
        return this;
    }

    @Override
    public Parameter getParameter() {
        return parameter;
    }

    @Override
    public int getDimension() {
        return parameter.getDimension();
    }

    @Override
    protected double calculateLogLikelihood() {
        double logLikelihood = 0.0;
        for (int i = 0; i < parameter.getDimension(); i++) {
            final double contribution = getLogDensityForCurrentCategory(i);
            if (!Double.isFinite(contribution)) {
                return Double.NEGATIVE_INFINITY;
            }
            logLikelihood += contribution;
        }
        return logLikelihood;
    }

    @Override
    public double[] getGradientLogDensity() {
        Arrays.fill(gradient, 0.0);
        for (int i = 0; i < parameter.getDimension(); i++) {
            final int nodeNumber = totalRewardsBranchRates.getNodeNumberFromParameterIndex(i);
            if (!rewardsAwareBranchModel.isAtomicBranch(nodeNumber)) {
                continue;
            }

            final int atomicState = rewardsAwareBranchModel.getAtomicBranchState(nodeNumber);
            gradient[i] = getGradientForAtomicState(i, atomicState);
        }
        return gradient;
    }

    public double getLogDensityForCurrentCategory(final int parameterIndex) {
        checkParameterIndex(parameterIndex);
        final int nodeNumber = totalRewardsBranchRates.getNodeNumberFromParameterIndex(parameterIndex);
        if (!rewardsAwareBranchModel.isAtomicBranch(nodeNumber)) {
            return getUniformLogDensity(parameterIndex);
        }
        return getLogDensityForAtomicState(
                parameterIndex,
                rewardsAwareBranchModel.getAtomicBranchState(nodeNumber));
    }

    public double getLogDensityForCategory(final int parameterIndex, final int category) {
        checkParameterIndex(parameterIndex);
        if (RewardMixtureCategoryDecoder.isContinuousCategory(category)) {
            return getUniformLogDensity(parameterIndex);
        }
        return getLogDensityForAtomicState(
                parameterIndex,
                RewardMixtureCategoryDecoder.getAtomicStateForCategory(category));
    }

    public double getStandardDeviation() {
        return standardDeviation;
    }

    @Override
    public double getNumericGradientStepSize() {
        return numericGradientStepSize;
    }

    @Override
    public void setNumericGradientStepSize(final double ratio) {
        numericGradientStepSize = ratio;
    }

    @Override
    public String getReport() {
        String report = "Reward-mixture atomic pseudo-prior; standardDeviation=" + standardDeviation + '\n';
        report += "analytic: " + Arrays.toString(getGradientLogDensity()) + '\n';
        return report;
    }

    private double getLogDensityForAtomicState(final int parameterIndex, final int atomicState) {
        final double x = parameter.getParameterValue(parameterIndex);
        final double mean = rewardsAwareBranchModel.getRewardRateRawForState(atomicState);
        final double lower = getLowerBound(parameterIndex);
        final double upper = getUpperBound(parameterIndex);

        if (!Double.isFinite(x) || x < lower || x > upper) {
            return Double.NEGATIVE_INFINITY;
        }

        final double logNormalDensity = NormalDistribution.logPdf(x, mean, standardDeviation);
        final double logTruncationMass = getLogTruncationMass(mean, lower, upper);
        if (!Double.isFinite(logTruncationMass)) {
            return Double.NEGATIVE_INFINITY;
        }
        return logNormalDensity - logTruncationMass;
    }

    private double getGradientForAtomicState(final int parameterIndex, final int atomicState) {
        final double x = parameter.getParameterValue(parameterIndex);
        final double mean = rewardsAwareBranchModel.getRewardRateRawForState(atomicState);
        return (mean - x) / variance;
    }

    private double getUniformLogDensity(final int parameterIndex) {
        final double x = parameter.getParameterValue(parameterIndex);
        final double lower = getLowerBound(parameterIndex);
        final double upper = getUpperBound(parameterIndex);
        if (!Double.isFinite(x) || x < lower || x > upper) {
            return Double.NEGATIVE_INFINITY;
        }
        return -Math.log(upper - lower);
    }

    private double getLogTruncationMass(final double mean, final double lower, final double upper) {
        final double logUpperCdf = NormalDistribution.cdf(upper, mean, standardDeviation, true);
        final double logLowerCdf = NormalDistribution.cdf(lower, mean, standardDeviation, true);
        return logSubtract(logUpperCdf, logLowerCdf);
    }

    private static double logSubtract(final double logA, final double logB) {
        if (Double.isInfinite(logB) && logB < 0.0) {
            return logA;
        }
        if (logB >= logA) {
            return Double.NEGATIVE_INFINITY;
        }
        return logA + Math.log1p(-Math.exp(logB - logA));
    }

    private double getLowerBound(final int parameterIndex) {
        return parameter.getBounds().getLowerLimit(parameterIndex);
    }

    private double getUpperBound(final int parameterIndex) {
        return parameter.getBounds().getUpperLimit(parameterIndex);
    }

    private void validateDimension() {
        final int branchCount = rewardsAwareBranchModel.getTree().getNodeCount() - 1;
        if (parameter.getDimension() != branchCount) {
            throw new IllegalArgumentException(
                    "totalRewardsBranchRates dimension must match branch count. Found " +
                            parameter.getDimension() + " but expected " + branchCount);
        }
    }

    private void validateBounds() {
        final Bounds<Double> bounds = parameter.getBounds();
        if (bounds.getBoundsDimension() != parameter.getDimension()) {
            throw new IllegalArgumentException(
                    "cts parameter bounds dimension must match parameter dimension. Found " +
                            bounds.getBoundsDimension() + " but expected " + parameter.getDimension());
        }
        for (int i = 0; i < parameter.getDimension(); i++) {
            final double lower = bounds.getLowerLimit(i);
            final double upper = bounds.getUpperLimit(i);
            if (!Double.isFinite(lower) || !Double.isFinite(upper) || !(upper > lower)) {
                throw new IllegalArgumentException(
                        "cts parameter bounds must be finite with upper > lower at index " + i +
                                ": lower=" + lower + ", upper=" + upper);
            }
        }
    }

    private void checkParameterIndex(final int parameterIndex) {
        if (parameterIndex < 0 || parameterIndex >= parameter.getDimension()) {
            throw new IllegalArgumentException("Parameter index out of range: " + parameterIndex);
        }
    }
}
